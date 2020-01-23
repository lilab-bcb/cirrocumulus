import logging
import os

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import scipy.sparse

logger = logging.getLogger("cirro")


def write_table(d, output_dir, name, write_statistics=True, row_group_size=None):
    os.makedirs(output_dir, exist_ok=True)
    pq.write_table(pa.Table.from_pydict(d), output_dir + os.path.sep + name + '.parquet',
        write_statistics=write_statistics, row_group_size=row_group_size)


def save_adata(adata, output_directory, X_range=None):
    logger.info('Save adata')
    if X_range is None:
        X_range = (0, adata.shape[1])
    logger.info('Save adata {}-{}'.format(X_range[0], X_range[1]))
    save_adata_X_chunk(adata, slice(X_range[0], X_range[1]), output_directory)
    if X_range[0] == 0:
        save_data_obsm(adata, output_directory)
        save_data_obs(adata, output_directory)


def save_data_obsm(adata, output_directory):
    logger.info('writing adata obsm')

    for name in adata.obsm.keys():
        m = adata.obsm[name]
        dimensions = [2]
        if m.shape[1] > 2:
            dimensions.append(3)
        d = {}
        for ndim in dimensions:
            for i in range(ndim):
                d[str(i)] = m[:, i].astype('float32')
            write_table(d, output_directory, name)


def save_data_obs(adata, output_directory):
    logger.info('writing adata obs')
    for name in adata.obs:
        # TODO sort?, row group size?
        write_table(dict(value=adata.obs[name]),
            output_directory, name)
    write_table(dict(value=adata.obs.index), output_directory, 'index')


def save_adata_X_chunk(adata, adata_col_slice, output_directory):
    X_slice = adata.X[:, adata_col_slice]
    names = adata.var.index[adata_col_slice]
    for j in range(X_slice.shape[1]):
        X = X_slice[:, j]
        if scipy.sparse.issparse(X):
            X = X.toarray().flatten()
        indices = np.where(X != 0)[0]
        values = X[indices]
        write_table(dict(index=indices, value=values), output_directory, names[j])
        if j > 0 and (j + 1) % 1000 == 0:
            logger.info('Wrote adata X {}/{}'.format(j + 1, X_slice.shape[1]))
