import gzip
import logging
import os

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import scipy.sparse

logger = logging.getLogger("cirro")


def write_pq(d, output_dir, name, write_statistics=True, row_group_size=None):
    os.makedirs(output_dir, exist_ok=True)
    pq.write_table(pa.Table.from_pydict(d), output_dir + os.path.sep + name + '.parquet',
        write_statistics=write_statistics, row_group_size=row_group_size)


def save_adata_pq(adata, schema, output_directory):
    import pandas._libs.json as ujson
    logger.info('Save adata')
    X_dir = os.path.join(output_directory, 'X')
    obs_dir = os.path.join(output_directory, 'obs')
    obsm_dir = os.path.join(output_directory, 'obsm')
    os.makedirs(X_dir, exist_ok=True)
    os.makedirs(obs_dir, exist_ok=True)
    os.makedirs(obsm_dir, exist_ok=True)
    with gzip.open(os.path.join(output_directory, 'index.json.gz'), 'wt') as f:
        # json.dump(result, f)
        f.write(ujson.dumps(schema, double_precision=2, orient='values'))

    save_adata_X(adata, X_dir)
    save_data_obs(adata, obs_dir)
    save_data_obsm(adata, obsm_dir)


def save_adata_X(adata, X_dir):
    adata_X = adata.X
    names = adata.var.index
    is_sparse = scipy.sparse.issparse(adata_X)
    if is_sparse and scipy.sparse.isspmatrix_csr(adata_X):
        adata_X = adata_X.tocsc()

    for j in range(adata_X.shape[1]):
        X = adata_X[:, j]
        if is_sparse:
            X = X.toarray().flatten()
        indices = np.where(X != 0)[0]
        values = X[indices]
        name = names[j]
        write_pq(dict(index=indices, value=values), X_dir, name)
        if j > 0 and (j + 1) % 1000 == 0:
            logger.info('Wrote adata X {}/{}'.format(j + 1, adata_X.shape[1]))


def save_data_obsm(adata, obsm_dir):
    logger.info('writing adata obsm')

    for name in adata.obsm.keys():
        m = adata.obsm[name]
        dimensions = 3 if m.shape[1] > 2 else 2
        d = {}
        for i in range(dimensions):
            d[name + '_' + str(i + 1)] = m[:, i].astype('float32')
        write_pq(d, obsm_dir, name)


def save_data_obs(adata, obs_dir):
    logger.info('writing adata obs')
    for name in adata.obs:
        # TODO sort?
        value = adata.obs[name]
        write_pq(dict(value=value), obs_dir, name)
    write_pq(dict(value=adata.obs.index.values), obs_dir, 'index')
