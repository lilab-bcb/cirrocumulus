import gzip

import numpy as np
import pandas as pd
import pandas._libs.json as ujson
import scipy.stats as ss

from .data_processing import get_filter_expr, data_filter_keys, get_type_to_measures
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


def get_mask(dataset_api, dataset, data_filter):
    measures, dimensions, basis_list = data_filter_keys(data_filter)
    keys = get_type_to_measures(measures)
    keys['obs'] = list(dimensions)
    keys['basis'] = basis_list
    df = dataset_api.read_dataset(keys=keys, dataset=dataset)
    return get_filter_expr(df, data_filter)


def run_job(database_api, dataset_api, email, job_id, job_type, dataset, params):
    database_api.update_job(email=email, job_id=job_id, status='running', result=None)
    schema = dataset_api.schema(dataset)
    var_names = schema['var']
    nfeatures = len(var_names)
    mask1 = get_mask(dataset_api, dataset, params['filter'])
    batch_size = 1000  # TODO
    scores = np.zeros(nfeatures)
    pvals = np.ones(nfeatures)
    is_sparse = None
    if job_type == 'de':
        fold_changes = np.zeros(nfeatures)
        mask2 = get_mask(dataset_api, dataset, params['filter2'])
        # overlapping cells are removed
    elif job_type == 'corr':
        df = dataset_api.read_dataset(keys=dict(X=[params['ref']]), dataset=dataset)
        df = df[mask1]
        ref = df[params['ref']]
        is_sparse = hasattr(ref, 'sparse')
        if is_sparse:
            ref = ref.sparse.to_dense()
        is_pearson = False  # params.get('method', 'pearson') == 'pearson'
        if not is_pearson:
            ref_expressed = (ref != 0)
            ref_not_expressed = ~ref_expressed
    index = 0
    for i in range(0, nfeatures, batch_size):
        start = i
        end = min(nfeatures, start + batch_size)
        features = var_names[start:end]
        df = dataset_api.read_dataset(keys=dict(X=features), dataset=dataset)
        df1 = df[mask1]
        if job_type == 'de':
            df2 = df[mask2]
            for feature in features:
                v1 = df1[feature]
                v2 = df2[feature]
                fc = v1.mean() - v2.mean()
                fold_changes[index] = fc
                if is_sparse is None:
                    is_sparse = hasattr(v1, 'sparse')
                if is_sparse:
                    v1 = v1.sparse.to_dense()
                    v2 = v2.sparse.to_dense()
                try:
                    scores[index], pvals[index] = ss.mannwhitneyu(v1, v2, alternative="two-sided")
                except ValueError:  # All numbers are identical
                    pass
                index += 1
        elif job_type == 'corr':
            for feature in features:
                v1 = df1[feature]
                if is_sparse:
                    v1 = v1.sparse.to_dense()
                try:
                    if is_pearson:
                        scores[index], pvals[index] = ss.pearsonr(ref, v1)
                    else:
                        #              gene1_exp, gene1_not_exp
                        # gene2_exp
                        # gene2_not_exp

                        a = (v1[ref_expressed] != 0).sum()
                        b = (v1[ref_not_expressed] != 0).sum()
                        c = (v1[ref_expressed] == 0).sum()
                        d = (v1[ref_not_expressed] == 0).sum()
                        scores[index], pvals[index] = ss.fisher_exact(
                            [[a, b], [c, d]])
                except ValueError:  # All numbers are identical
                    pass
                index += 1
        database_api.update_job(email=email, job_id=job_id, status='running {:.0f}%'.format(100 * index / nfeatures),
            result=None)
    pvals = fdrcorrection(pvals)
    if job_type == 'de':
        result_df = pd.DataFrame(
            data={'index': var_names, 'comparison:pvals_adj': pvals, 'comparison:fold_changes': fold_changes,
                  'comparison:scores': scores})
        result = dict(groups=['comparison'], fields=['pvals_adj', 'fold_changes', 'scores'],
            data=result_df.to_dict(orient='records'))
    elif job_type == 'corr':
        result_df = pd.DataFrame(
            data={'index': var_names, 'selection:pvals_adj': pvals, 'selection:scores': scores})
        result = dict(groups=['selection'], fields=['pvals_adj', 'scores'], data=result_df.to_dict(orient='records'))
    result = ujson.dumps(result, double_precision=2, orient='values').encode('UTF-8')
    result = gzip.compress(result)
    database_api.update_job(email=email, job_id=job_id, status='complete', result=result)
