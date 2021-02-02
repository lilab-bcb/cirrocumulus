import gzip
from concurrent.futures.thread import ThreadPoolExecutor

import pandas as pd
import pandas._libs.json as ujson

from .diff_exp import fdrcorrection

executor = ThreadPoolExecutor(max_workers=2)


def submit_job(database_api, dataset_api, email, dataset, job_name, job_type, params):
    job_id = database_api.create_job(email=email, dataset_id=dataset['id'], job_name=job_name, job_type=job_type,
        params=params)
    f = diff_exp
    # executor.submit(f, database_api, dataset_api, email, job_id, dataset, params)
    f(database_api, dataset_api, email, job_id, dataset, params)
    return job_id


def diff_exp(database_api, dataset_api, email, job_id, dataset, params):
    database_api.update_job(email=email, job_id=job_id, status='running', result=None)
    import numpy as np
    import scipy.stats as ss
    schema = dataset_api.schema(dataset)
    var_names = schema['var']
    shape = schema['shape']
    nfeatures = len(var_names)
    print(nfeatures, shape)
    scores = np.full(nfeatures, 0)
    pvals = np.full(nfeatures, 1)
    fold_changes = np.full(nfeatures, 0)
    batch_size = 100
    index = 0
    mask1 = np.random.choice(a=[False, True], size=schema['shape'][0], p=[0.5, 0.5])  # FIXME
    mask2 = ~mask1

    for i in range(0, nfeatures, batch_size):
        start = i
        end = min(nfeatures, start + batch_size)
        features = var_names[start:end]
        df = dataset_api.read_dataset(keys=dict(X=features), dataset=dataset)
        df1 = df[mask1]
        df2 = df[mask2]
        for feature in features:
            v1 = df1[feature]
            v2 = df2[feature]
            fc = v1.mean() - v2.mean()
            fold_changes[index] = fc
            is_sparse = hasattr(v1, 'sparse')
            if is_sparse:
                v1 = v1.sparse.to_dense()
                v2 = v2.sparse.to_dense()
            if is_sparse:
                try:
                    scores[index], pvals[index] = ss.mannwhitneyu(v1.values, v2.values,
                        alternative="two-sided")
                except ValueError:
                    # All numbers are identical
                    pass
            else:
                try:
                    scores[index], pvals[index] = ss.mannwhitneyu(v1, v2, alternative="two-sided")
                except ValueError:
                    # All numbers are identical
                    pass
            index += 1

    pvals = fdrcorrection(pvals)
    fields = ['pvals_adj', 'fold_changes', 'scores']
    result_df = pd.DataFrame(
        data={'index': var_names, '1:pvals_adj': pvals, '1:fold_changes': fold_changes, '1:scores': scores})
    result = dict(groups=['1'], fields=fields, data=result_df.to_dict(orient='records'))
    result = ujson.dumps(result, double_precision=2, orient='values').encode('UTF-8')
    result = gzip.compress(result)
    database_api.update_job(email=email, job_id=job_id, status='complete', result=result)
