import os
import math
import logging

import pandas as pd

from cirrocumulus.diff_exp import DE
from cirrocumulus.ot.transport_map_model import read_transport_map_dir
from cirrocumulus.util import dumps

from .data_processing import get_filter_str, get_mask, get_selected_data
from .envir import (
    CIRRO_DATABASE_CLASS,
    CIRRO_JOB_RESULTS,
    CIRRO_JOB_TYPE,
    CIRRO_MAX_WORKERS,
    CIRRO_SERVE,
)
from .fdr import fdrcorrection
from .util import add_dataset_providers, create_instance, get_fs, import_path, open_file


executor = None
job_id_2_future = dict()

logger = logging.getLogger("cirro")


def save_job_result_to_file(result, job_id):
    new_result = dict()
    new_result["content-type"] = result.pop("content-type")
    if new_result["content-type"] == "application/json":
        new_result["content-encoding"] = "gzip"
        url = os.path.join(os.environ[CIRRO_JOB_RESULTS], str(job_id) + ".json.gz")
        with open_file(url, "wt", compression="gzip") as out:
            out.write(dumps(result, double_precision=2, orient="values"))
    elif new_result["content-type"] == "application/h5ad":
        url = os.path.join(os.environ[CIRRO_JOB_RESULTS], str(job_id) + ".h5ad")
        with get_fs(url).open(url, "wb") as out:
            result["data"].write(out)
    elif new_result["content-type"] == "application/zarr":
        url = os.path.join(os.environ[CIRRO_JOB_RESULTS], str(job_id) + ".zarr")
        result["data"].write_zarr(get_fs(url).get_mapper(url))
    elif new_result["content-type"] == "application/parquet":
        import pyarrow as pa
        import pyarrow.parquet as pq

        url = os.path.join(os.environ[CIRRO_JOB_RESULTS], str(job_id) + ".parquet")
        pq.write_table(pa.Table.from_pandas(result["data"]), url, filesystem=get_fs(url))
    else:
        raise ValueError("Unknown content-type {}".format(new_result["content-type"]))
    new_result["url"] = url
    return new_result


def delete_job(job_id):
    future = job_id_2_future.get(job_id)
    if future is not None and not future.done():
        del job_id_2_future[job_id]
        future.cancel()
        logger.info("Cancel job {}".format(job_id))


def done_callback(future):
    for job_id in list(job_id_2_future.keys()):
        if job_id_2_future[job_id] == future:
            del job_id_2_future[job_id]
            logger.info("Job {} done".format(job_id))
            break


def submit_job(database_api, dataset_api, email, dataset, job_name, job_type, params):
    global executor

    is_serve = os.environ.get(CIRRO_SERVE) == "true"
    if executor is None:
        max_workers = int(os.environ.get(CIRRO_MAX_WORKERS, "2" if is_serve else "1"))
        if max_workers > 0:
            import multiprocessing
            from concurrent.futures.process import ProcessPoolExecutor
            from concurrent.futures.thread import ThreadPoolExecutor

            executor = (
                ProcessPoolExecutor(
                    max_workers=max_workers, mp_context=multiprocessing.get_context("spawn")
                )
                if is_serve
                else ThreadPoolExecutor(max_workers=max_workers)
            )
    job_id = database_api.create_job(
        email=email, dataset_id=dataset["id"], job_name=job_name, job_type=job_type, params=params
    )
    if executor is not None:
        future = executor.submit(
            run_job,
            email,
            job_id,
            job_name,
            job_type,
            dataset,
            params,
            database_api if not is_serve else None,
            dataset_api if not is_serve else None,
        )
        future.add_done_callback(done_callback)
        job_id_2_future[job_id] = future
    else:
        run_job(
            email,
            job_id,
            job_name,
            job_type,
            dataset,
            params,
            database_api if not is_serve else None,
            dataset_api if not is_serve else None,
        )
    return job_id


def get_obs(dataset_api, dataset, dataset_info, params):
    obs_fields = params.get("obs")
    if obs_fields is not None:
        obs = dataset_api.read_dataset(keys=dict(obs=obs_fields), dataset=dataset).obs
        obs_field = obs_fields[0]
        if len(obs_fields) > 1:
            # combine in to one field
            obs_field = "_".join(obs_fields)
            obs[obs_field] = obs[obs_fields[0]].astype(str)
            for i in range(1, len(obs_fields)):
                obs[obs_field] += "_" + obs[obs_fields[i]].astype(str)
            obs[obs_field] = obs[obs_field].astype("category")
        return obs, obs_field
    else:
        filters = [params["filter"], params["filter2"]]
        filter_names = [get_filter_str(params["filter"]), get_filter_str(params["filter2"])]

        for i in range(len(filter_names)):
            if filter_names[i] is None:
                filter_names[i] = "group_" + str(i + 1)
        if filter_names[1] == filter_names[0]:
            filter_names[1] = f"{filter_names[1]}_2"
        obs = pd.DataFrame(index=pd.RangeIndex(dataset_info["shape"][0]).astype(str))
        obs_field = "selection"  # order of categories needs to match filter names
        # obs[obs_field] = "3"
        masks, _ = get_mask(dataset_api, dataset, dataset_info, filters)
        for i in range(2):
            obs.loc[masks[i], obs_field] = filter_names[i]
        obs[obs_field] = obs[obs_field].astype(pd.Categorical(filter_names, ordered=True).dtype)
        return obs, obs_field


def run_job(email, job_id, job_name, job_type, dataset, params, database_api, dataset_api):
    if database_api is None:
        database_api = create_instance(os.environ[CIRRO_DATABASE_CLASS])
    if dataset_api is None:
        from cirrocumulus.api import dataset_api

        add_dataset_providers()
    database_api.update_job(email=email, job_id=job_id, status="running", result=None)
    f = os.environ[CIRRO_JOB_TYPE + job_type]
    if f is None:
        database_api.update_job(email=email, job_id=job_id, status="error", result=None)
        raise ValueError("No function to handle {}".format(job_type))
    import_path(f)(email, job_id, job_name, job_type, dataset, params, database_api, dataset_api)


def run_ot_trajectory(
    email, job_id, job_name, job_type, dataset, params, database_api, dataset_api
):
    dataset_info = dataset_api.get_dataset_info(dataset)
    timepoint_field = dataset_info.get("timepoint_field", "day")
    selected_adata = get_selected_data(
        dataset_api,
        dataset,
        measures=["obs/index", "obs/" + timepoint_field],
        data_filter=params["filter"],
        dataset_info=dataset_info,
    )

    day = selected_adata.obs[timepoint_field].unique()
    if len(day) > 1:
        database_api.update_job(
            email=email, job_id=job_id, status="error", result="More than one timepoint selected"
        )
        raise ValueError("More than one timepoint selected")
    tmap_name = params["tmap"]
    for tmap_item in dataset["ot"]["tmaps"]:
        if tmap_item["name"] == tmap_name:
            break

    tmap_model = read_transport_map_dir(tmap_item["path"])
    populations = tmap_model.population_from_ids(selected_adata.obs["index"].values, at_time=day[0])
    populations[0].name = job_name

    trajectory_dataset = tmap_model.trajectories(populations)

    # match order in dataset
    # trajectory_ds = trajectory_ds[adata.obs.index]

    result = dict(data=trajectory_dataset.to_df())
    result["content-type"] = "application/parquet"
    database_api.update_job(email=email, job_id=job_id, status="complete", result=result)


def run_de(email, job_id, job_name, job_type, dataset, params, database_api, dataset_api):
    dataset_info = dataset_api.get_dataset_info(dataset)
    var_names = dataset_info["var"]
    nfeatures = len(var_names)
    batch_size = (
        math.ceil(nfeatures / 5) if os.environ.get(CIRRO_SERVE) == "true" else nfeatures
    )  # TODO more intelligent batches

    def get_batch_fn(i):
        end = min(nfeatures, i + batch_size)
        adata = dataset_api.read_dataset(keys=dict(X=[slice(i, end)]), dataset=dataset)
        if batch_size != nfeatures:
            frac = end / nfeatures
            status = "running {:.0f}%".format(100 * frac) if frac < 1 else "saving results"
            logger.info(status)
            database_api.update_job(email=email, job_id=job_id, status=status, result=None)
        return adata

    compare_pairs = params.get("pairs") == "1"
    obs, obs_field = get_obs(dataset_api, dataset, dataset_info, params)

    de_results = DE(
        series=obs[obs_field],
        nfeatures=nfeatures,
        batch_size=batch_size,
        get_batch_fn=get_batch_fn,
        one_vs_rest=not compare_pairs,
    )  # TODO get base

    # group:field is object entry
    result_df = pd.DataFrame(data={"index": var_names})
    has_frac_expressed = False
    comparison_names = []
    for comparison in de_results.pair2results.keys():
        result = de_results.pair2results[comparison]
        comparison_name = comparison if isinstance(comparison, str) else "_".join(comparison)
        comparison_names.append(comparison_name)

        pvals = fdrcorrection(result["pvals"])
        result_df[f"{comparison_name}:pvals_adj"] = pvals
        result_df[f"{comparison_name}:scores"] = result["scores"]
        result_df[f"{comparison_name}:lfc"] = result["logfoldchanges"]

        if result.get("frac_expressed1") is not None:
            has_frac_expressed = True
            result_df[f"{comparison_name}:pts_1"] = result["frac_expressed1"]
            result_df[f"{comparison_name}:pts_2"] = result["frac_expressed2"]
    # client expects field {comparison_name}:pvals_adj
    result = dict(
        groups=comparison_names,
        fields=["pvals_adj", "scores", "lfc"] + (["pts_1", "pts_2"] if has_frac_expressed else []),
        data=result_df.to_dict(orient="records"),
    )
    result["content-type"] = "application/json"
    database_api.update_job(email=email, job_id=job_id, status="complete", result=result)
