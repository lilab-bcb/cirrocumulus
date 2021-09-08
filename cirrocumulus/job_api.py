import numpy as np
import pandas as pd
from anndata import AnnData
from cirrocumulus.de import DE

from .data_processing import get_filter_str, get_mask
from .diff_exp import fdrcorrection
from .envir import CIRRO_SERVE, CIRRO_MAX_WORKERS

executor = None


def submit_job(database_api, dataset_api, email, dataset, job_name, job_type, params):
    global executor
    if executor is None:
        from concurrent.futures.thread import ThreadPoolExecutor
        import os

        if os.environ.get(CIRRO_SERVE) != 'true':
            max_workers = 1
        else:
            max_workers = int(os.environ.get(CIRRO_MAX_WORKERS, '2'))
        executor = ThreadPoolExecutor(max_workers=max_workers)
    job_id = database_api.create_job(email=email, dataset_id=dataset['id'], job_name=job_name, job_type=job_type,
                                     params=params)
    # run_job(database_api, dataset_api, email, job_id, job_type, dataset, params)
    executor.submit(run_job, database_api, dataset_api, email, job_id, job_type, dataset, params)
    return job_id


def _power(X, power):
    return X ** power if isinstance(X, np.ndarray) else X.power(power)


def run_job(database_api, dataset_api, email, job_id, job_type, dataset, params):
    database_api.update_job(email=email, job_id=job_id, status='running', result=None)
    schema = dataset_api.schema(dataset)
    var_names = schema['var']
    nfeatures = len(var_names)
    filters = [params['filter']]
    if job_type == 'de':
        filters.append(params['filter2'])
    masks, _ = get_mask(dataset_api, dataset, filters)
    batch_size = 5000  # TODO more intelligent batches
    obs_field = 'tmp'

    def get_batch_fn(i):
        start = i
        end = min(nfeatures, start + batch_size)
        features = var_names[start:end]
        adata = dataset_api.read_dataset(keys=dict(X=features), dataset=dataset)
        if batch_size != nfeatures:
            database_api.update_job(email=email, job_id=job_id, status='running {:.0f}%'.format(100 * end / nfeatures),
                                    result=None)
        return adata

    obs = pd.DataFrame(index=pd.RangeIndex(schema['shape'][0]))
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
                  format='json',
                  fields=['pvals_adj', 'scores', 'lfc'] + (['pts_1', 'pts_2'] if de.frac_expressed is not None else []),
                  data=result_df.to_dict(orient='records'))
    database_api.update_job(email=email, job_id=job_id, status='complete', result=result)
