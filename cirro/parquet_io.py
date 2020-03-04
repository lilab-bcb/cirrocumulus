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


def save_adata(adata, root, X_range=None):
    logger.info('Save adata')
    if X_range is None:
        X_range = (0, adata.shape[1])
    logger.info('Save adata {}-{}'.format(X_range[0], X_range[1]))
    X_dir = os.makedirs(os.path.join(root, 'X'), exist_ok=True)
    os.makedirs(os.path.join(root, 'obs'))
    os.makedirs(os.path.join(root, 'obsm'))

    save_adata_X_chunk(adata, slice(X_range[0], X_range[1]), root)
    if X_range[0] == 0:
        import json
        import os
        import gzip
        with gzip.open(os.path.join(X_dir, 'index.json.gz'), 'wt') as f:
            json.dump({'shape': adata.shape, 'encoding': 'parquet', 'version': '1.0.0'}, f)
        save_data_obs(adata, root)
        save_data_obsm(adata, root)


def save_adata_X_chunk(adata, adata_col_slice, X_dir):
    X_slice = adata.X[:, adata_col_slice]
    names = adata.var.index[adata_col_slice]
    is_sparse = scipy.sparse.issparse(X_slice)
    if is_sparse and scipy.sparse.isspmatrix_csr(X_slice):
        X_slice = X_slice.tocsc()

    for j in range(X_slice.shape[1]):
        X = X_slice[:, j]
        if is_sparse:
            X = X.toarray().flatten()
        indices = np.where(X != 0)[0]
        values = X[indices]
        write_pq(dict(index=indices, value=values), os.path.join(X_dir, 'name=' + names[j]), 'index')
        if j > 0 and (j + 1) % 1000 == 0:
            logger.info('Wrote adata X {}/{}'.format(j + 1, X_slice.shape[1]))


def save_data_obsm(adata, obsm_dir):
    logger.info('writing adata obsm')

    for name in adata.obsm.keys():
        m = adata.obsm[name]
        dimensions = 3 if m.shape[1] > 2 else 2
        d = {}
        for i in range(dimensions):
            d[name + '_' + str(i + 1)] = m[:, i].astype('float32')
        write_pq(d, os.path.join(obsm_dir, name), 'index.parquet')


def save_data_obs(adata, obs_dir):
    logger.info('writing adata obs')
    for name in adata.obs:
        # TODO sort?
        value = adata.obs[name]
        write_pq(dict(value=value), os.path.join(obs_dir, name), 'index')
    write_pq(dict(value=adata.obs.index.values), output, 'index')
