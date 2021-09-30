import logging
import os

import numpy as np
import pandas as pd
import pandas._libs.json as ujson
from anndata import AnnData

from cirrocumulus.de import DE
from .data_processing import get_filter_str, get_mask
from .diff_exp import fdrcorrection
from .envir import CIRRO_SERVE, CIRRO_MAX_WORKERS, CIRRO_DATABASE_CLASS, CIRRO_JOB_RESULTS
from .util import create_instance, add_dataset_providers, get_fs

executor = None

logger = logging.getLogger('cirro')
job_type_to_func = dict()


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


def submit_job(database_api, dataset_api, email, dataset, job_name, job_type, params):
    global executor
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
    executor.submit(run_job, email, job_id, job_type, dataset, params, database_api if not is_serve else None,
                    dataset_api if not is_serve else None)
    return job_id


def _power(X, power):
    return X ** power if isinstance(X, np.ndarray) else X.power(power)


def run_job(email, job_id, job_type, dataset, params, database_api, dataset_api):
    if database_api is None:
        database_api = create_instance(os.environ[CIRRO_DATABASE_CLASS])
    if dataset_api is None:
        from cirrocumulus.api import dataset_api
        add_dataset_providers()
    database_api.update_job(email=email, job_id=job_id, status='running', result=None)
    f = job_type_to_func[job_type]
    f(email, job_id, job_type, dataset, params, database_api, dataset_api)


def run_de(email, job_id, job_type, dataset, params, database_api, dataset_api):
    dataset_info = dataset_api.get_dataset_info(dataset)
    var_names = dataset_info['var']
    nfeatures = len(var_names)
    filters = [params['filter']]
    if job_type == 'de':
        filters.append(params['filter2'])
    masks, _ = get_mask(dataset_api, dataset, filters)
    batch_size = 5000 if os.environ.get(CIRRO_SERVE) == 'true' else nfeatures  # TODO more intelligent batches
    obs_field = 'tmp'

    def get_batch_fn(i):
        logger.info('batch {}'.format(i))
        end = min(nfeatures, i + batch_size)
        adata = dataset_api.read_dataset(keys=dict(X=[slice(i, end)]), dataset=dataset)
        if batch_size != nfeatures:
            frac = end / nfeatures
            database_api.update_job(email=email, job_id=job_id,
                                    status='running {:.0f}%'.format(100 * frac) if frac < 1 else 'saving results',
                                    result=None)
        return adata

    obs = pd.DataFrame(index=pd.RangeIndex(dataset_info['shape'][0]))
    obs[obs_field] = ''
    for i in range(len(masks)):
        obs.loc[masks[i], obs_field] = str(i)
    obs[obs_field] = obs[obs_field].astype('category')
    de = DE(AnnData(obs=obs), obs_field, nfeatures, batch_size, get_batch_fn,
            key_set=[str(i) for i in range(len(masks))])
    pvals = fdrcorrection(de.pvals)
    # group:field is object entry
    if job_type == 'de':
        comparison1_name = get_filter_str(params['filter'])
    comparison2_name = get_filter_str(params['filter2'])
    comparison = 'comparison' if comparison1_name is None or comparison2_name is None else comparison1_name + '_' + comparison2_name

    result_df = pd.DataFrame(
        data={'index': var_names, '{}:pvals_adj'.format(comparison): pvals, '{}:scores'.format(comparison): de.scores,
              '{}:lfc'.format(comparison): de.logfoldchanges})
    if de.frac_expressed is not None:
        result_df['{}:pts_1'.format(comparison)] = de.frac_expressed.iloc[0]
        result_df['{}:pts_2'.format(comparison)] = de.frac_expressed.iloc[1]
    result = dict(groups=[comparison],
                  fields=['pvals_adj', 'scores', 'lfc'] + (['pts_1', 'pts_2'] if de.frac_expressed is not None else []),
                  data=result_df.to_dict(orient='records'))
    result['content-type'] = 'application/json'
    database_api.update_job(email=email, job_id=job_id, status='complete', result=result)


job_type_to_func['de'] = run_de
