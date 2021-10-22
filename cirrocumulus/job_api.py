import logging
import os

import anndata
import math
import numpy as np
import pandas as pd
import pandas._libs.json as ujson
from anndata import AnnData

from cirrocumulus.de import DE
from .data_processing import get_filter_str, get_mask
from .diff_exp import fdrcorrection
from .envir import CIRRO_SERVE, CIRRO_MAX_WORKERS, CIRRO_DATABASE_CLASS, CIRRO_JOB_RESULTS, CIRRO_JOB_TYPE
from .util import create_instance, add_dataset_providers, get_fs, import_path

executor = None
job_id_2_future = dict()

logger = logging.getLogger('cirro')


def save_job_result_to_file(result, job_id):
    new_result = dict()
    new_result['content-type'] = result.pop('content-type')
    if new_result['content-type'] == 'application/json':
        compression = 'gzip'
        new_result['content-encoding'] = compression
        url = os.path.join(os.environ[CIRRO_JOB_RESULTS], str(job_id) + '.json.gz')
        with get_fs(url).open(url, 'wt', compression=compression) as out:
            out.write(ujson.dumps(result, double_precision=2, orient='values'))
    elif new_result['content-type'] == 'application/h5ad':
        url = os.path.join(os.environ[CIRRO_JOB_RESULTS], str(job_id) + '.h5ad')
        with get_fs(url).open(url, 'wb') as out:
            result['data'].write(out)
    elif new_result['content-type'] == 'application/zarr':
        url = os.path.join(os.environ[CIRRO_JOB_RESULTS], str(job_id) + '.zarr')
        result['data'].write_zarr(get_fs(url).get_mapper(url))
    elif new_result['content-type'] == 'application/parquet':
        import pyarrow.parquet as pq
        import pyarrow as pa
        url = os.path.join(os.environ[CIRRO_JOB_RESULTS], str(job_id) + '.parquet')
        pq.write_table(pa.Table.from_pandas(result['data']), url, filesystem=get_fs(url))
    else:
        raise ValueError('Unknown content-type {}'.format(new_result['content-type']))
    new_result['url'] = url
    return new_result


def delete_job(job_id):
    global job_id_2_future
    future = job_id_2_future.get(job_id)
    if future is not None and not future.done():
        del job_id_2_future[job_id]
        future.cancel()
        logger.info('Cancel job {}'.format(job_id))


def done_callback(future):
    global job_id_2_future
    for job_id in list(job_id_2_future.keys()):
        if job_id_2_future[job_id] == future:
            del job_id_2_future[job_id]
            logger.info('Job {} done'.format(job_id))
            break


def submit_job(database_api, dataset_api, email, dataset, job_name, job_type, params):
    global executor
    global job_id_2_future
    from concurrent.futures.process import ProcessPoolExecutor
    from concurrent.futures.thread import ThreadPoolExecutor
    import os
    is_serve = os.environ.get(CIRRO_SERVE) == 'true'
    if executor is None:

        if not is_serve:
            max_workers = 1
        else:
            max_workers = int(os.environ.get(CIRRO_MAX_WORKERS, '2'))
        executor = ProcessPoolExecutor(max_workers=max_workers) if is_serve else ThreadPoolExecutor(
            max_workers=max_workers)
    job_id = database_api.create_job(email=email, dataset_id=dataset['id'], job_name=job_name, job_type=job_type,
                                     params=params)
    # run_job(email, job_id, job_type, dataset, params, database_api if not is_serve else None,
    #         dataset_api if not is_serve else None)
    future = executor.submit(run_job, email, job_id, job_type, dataset, params, database_api if not is_serve else None,
                             dataset_api if not is_serve else None)
    future.add_done_callback(done_callback)
    job_id_2_future[job_id] = future

    return job_id


def _power(X, power):
    return X ** power if isinstance(X, np.ndarray) else X.power(power)


def get_comparisons(dataset_api, dataset, dataset_info, params, X, combinations):
    obs_fields = params.get('obs')
    if obs_fields is not None:
        import itertools
        adata = dataset_api.read_dataset(keys=dict(obs=obs_fields, X=X if X is not None else []), dataset=dataset)
        obs_field = obs_fields[0]
        if len(obs_fields) > 1:
            # combine in to one field
            obs_field = '_'.join(obs_fields)
            adata.obs[obs_field] = adata.obs[obs_fields[0]].astype(str)
            for i in range(1, len(obs_fields)):
                adata.obs[obs_field] += '_' + adata.obs[obs_fields[i]].astype(str)
            adata.obs[obs_field] = adata.obs[obs_field].astype('category')
        comparison_names = list(
            itertools.combinations(adata.obs[obs_field].cat.categories, 2)) if combinations else list(
            itertools.permutations(adata.obs[obs_field].cat.categories, 2))
        return dict(adata=adata, comparison_names=comparison_names, obs_field=obs_field, is_single=False)
    else:
        filters = [params['filter'], params['filter2']]
        filter_names = [get_filter_str(params['filter']), get_filter_str(params['filter2'])]
        for i in range(len(filter_names)):
            if filter_names[i] is None:
                filter_names[i] = 'group_' + str(i + 1)
        obs = pd.DataFrame(index=pd.RangeIndex(dataset_info['shape'][0]))
        obs_field = 'user'
        obs[obs_field] = ''
        masks, _ = get_mask(dataset_api, dataset, filters)
        for i in range(len(masks)):
            obs.loc[masks[i], obs_field] = filter_names[i]
        obs[obs_field] = obs[obs_field].astype('category')
        if X is not None and len(X) > 0:
            adata = dataset_api.read_dataset(keys=dict(X=X), dataset=dataset)
            adata.obs = obs
        else:
            adata = anndata.AnnData(obs=obs)

        return dict(adata=adata, comparison_names=[tuple(filter_names)], obs_field=obs_field,
                    is_single=True)


def run_job(email, job_id, job_type, dataset, params, database_api, dataset_api):
    if database_api is None:
        database_api = create_instance(os.environ[CIRRO_DATABASE_CLASS])
    if dataset_api is None:
        from cirrocumulus.api import dataset_api
        add_dataset_providers()
    database_api.update_job(email=email, job_id=job_id, status='running', result=None)
    f = os.environ[CIRRO_JOB_TYPE + job_type]
    if f is None:
        database_api.update_job(email=email, job_id=job_id, status='error', result=None)
        raise ValueError('No function to handle {}'.format(job_type))
    import_path(f)(email, job_id, job_type, dataset, params, database_api, dataset_api)


def run_de(email, job_id, job_type, dataset, params, database_api, dataset_api):
    dataset_info = dataset_api.get_dataset_info(dataset)
    var_names = dataset_info['var']
    nfeatures = len(var_names)

    batch_size = math.ceil(nfeatures / 3) if os.environ.get(
        CIRRO_SERVE) == 'true' else nfeatures  # TODO more intelligent batches

    def get_batch_fn(i):
        end = min(nfeatures, i + batch_size)
        adata = dataset_api.read_dataset(keys=dict(X=[slice(i, end)]), dataset=dataset)
        if batch_size != nfeatures:
            frac = end / nfeatures
            status = 'running {:.0f}%'.format(100 * frac) if frac < 1 else 'saving results'
            logger.info(status)
            database_api.update_job(email=email, job_id=job_id, status=status, result=None)
        return adata

    comparison_result = get_comparisons(dataset_api, dataset, dataset_info, params, X=None, combinations=True)
    comparison_names = comparison_result['comparison_names']
    de = DE(AnnData(obs=comparison_result['adata'].obs), comparison_result['obs_field'], nfeatures, batch_size,
            get_batch_fn, key_set=comparison_names[0] if comparison_result['is_single'] else None,
            pairs=comparison_names)
    # group:field is object entry
    result_df = pd.DataFrame(data={'index': var_names})
    has_frac_expressed = False
    comparison_name_strs = []
    for comparison_name in comparison_names:
        result = de.pair2results[comparison_name]
        comparison_name = '_'.join(comparison_name)
        comparison_name_strs.append(comparison_name)
        pvals = fdrcorrection(result['pvals'])
        result_df[f'{comparison_name}:pvals_adj'] = pvals
        result_df[f'{comparison_name}:scores'] = result['scores']
        result_df[f'{comparison_name}:lfc'] = result['logfoldchanges']

        if result.get('frac_expressed1') is not None:
            has_frac_expressed = True
            result_df[f'{comparison_name}:pts_1'] = result['frac_expressed1']
            result_df[f'{comparison_name}:pts_2'] = result['frac_expressed2']
    # client expects field {comparison_name}:pvals_adj
    result = dict(groups=comparison_name_strs,
                  fields=['pvals_adj', 'scores', 'lfc'] + (
                      ['pts_1', 'pts_2'] if has_frac_expressed else []),
                  data=result_df.to_dict(orient='records'))
    result['content-type'] = 'application/json'
    database_api.update_job(email=email, job_id=job_id, status='complete', result=result)
