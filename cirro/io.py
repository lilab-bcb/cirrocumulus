import logging

import numpy as np
import scipy.sparse

from cirro.parquet_dataset import write_pq


logger = logging.getLogger("cirro")


def save_adata(adata, root, X_range=None):
    logger.info('Save adata')
    if X_range is None:
        X_range = (0, adata.shape[1])
    logger.info('Save adata {}-{}'.format(X_range[0], X_range[1]))
    save_adata_X_chunk(adata, slice(X_range[0], X_range[1]), root)
    if X_range[0] == 0:
        import json
        import os
        import gzip
        with gzip.open(os.path.join(root, 'index.json.gz'), 'wt') as f:
            json.dump({'shape': adata.shape, 'encoding': 'parquet', 'version': '1.0.0'}, f)
        save_data_obs(adata, root)
        save_data_obsm(adata, root)


def save_data_obsm(adata, output):
    logger.info('writing adata obsm')

    for name in adata.obsm.keys():
        m = adata.obsm[name]
        dimensions = 3 if m.shape[1] > 2 else 2
        d = {}
        for i in range(dimensions):
            d[name + '_' + str(i + 1)] = m[:, i].astype('float32')
        write_pq(d, output, name)


def save_data_obs(adata, output):
    logger.info('writing adata obs')
    for name in adata.obs:
        # TODO sort?
        value = adata.obs[name]
        write_pq(dict(value=value), output, name)
    write_pq(dict(value=adata.obs.index.values), output, 'index')


def save_adata_X_chunk(adata, adata_col_slice, output):
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
        write_pq(dict(index=indices, value=values), output, names[j])
        if j > 0 and (j + 1) % 1000 == 0:
            logger.info('Wrote adata X {}/{}'.format(j + 1, X_slice.shape[1]))
