from concurrent.futures.thread import ThreadPoolExecutor

from .diff_exp import fdrcorrection

executor = ThreadPoolExecutor(max_workers=2)


def submit_job(database_api, dataset_api, email, dataset, params):
    job_id = database_api.create_job(email=email, dataset_id=dataset['id'], params=params)
    executor.submit(diff_exp, database_api, dataset_api, email, job_id, dataset, params)
    return job_id


def diff_exp(database_api, dataset_api, email, job_id, dataset, params):
    database_api.update_job(email=email, job_id=job_id, status='running', result=None)
    import numpy as np
    import scipy.stats as ss
    schema = dataset_api.schema(dataset)
    var_names = schema['var']
    nfeatures = len(var_names)
    pvals = np.full(nfeatures, 1.0)
    fold_changes = np.full(nfeatures, 0)
    batch_size = 100
    index = 0
    mask1 = np.random.choice(a=[False, True], size=schema['shape'][0], p=[0.5, 0.5])  # FIXME
    mask2 = ~mask1
    import time
    start = time.time()
    for i in range(0, nfeatures, batch_size):
        start = i
        end = min(nfeatures, start + batch_size)
        features = var_names[start:end]
        df = dataset_api.read_dataset(keys=dict(X=features), dataset=dataset)
        for feature in features:
            v1 = df[mask1][feature]
            v2 = df[mask2][feature]
            fc = v1.mean() - v2.mean()
            fold_changes[index] = fc
            is_sparse = hasattr(v1, 'sparse')
            if is_sparse:
                v1 = v1.sparse.to_dense()
                v2 = v2.sparse.to_dense()
            if is_sparse:
                try:
                    _, pvals[index] = ss.mannwhitneyu(v1.values, v2.values,
                        alternative="two-sided")
                except ValueError:
                    # All numbers are identical
                    pass
            else:
                try:
                    _, pvals[index] = ss.mannwhitneyu(v1, v2, alternative="two-sided")
                except ValueError:
                    # All numbers are identical
                    pass
            index += 1

    print(time.time() - start)
    qvals = fdrcorrection(pvals)

    result = dict(qvals=qvals.tolist(), pvals=pvals.tolist(), fold_changes=fold_changes.tolist())
    database_api.update_job(email=email, job_id=job_id, status='complete', result=result)
    print('updated')
