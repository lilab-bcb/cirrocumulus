import logging
import os

import numcodecs
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import scipy.sparse
import zarr

logger = logging.getLogger("cirro")


def write_table(d, output_dir, name, write_statistics=True, row_group_size=None):
    os.makedirs(output_dir, exist_ok=True)
    pq.write_table(pa.Table.from_pydict(d), output_dir + os.path.sep + name + '.parquet',
        write_statistics=write_statistics, row_group_size=row_group_size)


def save_adata(adata, output_directory, X_range=None, output_format='parquet'):
    logger.info('Save adata')
    if X_range is None:
        X_range = (0, adata.shape[1])
    logger.info('Save adata {}-{}'.format(X_range[0], X_range[1]))
    if output_format == 'zarr':
        output_directory = zarr.open(output_directory, mode='w')
    save_adata_X_chunk(adata, slice(X_range[0], X_range[1]), output_directory)
    if X_range[0] == 0:
        save_data_obs(adata, output_directory)
        save_data_obsm(adata, output_directory)


def save_data_obsm(adata, output):
    logger.info('writing adata obsm')
    is_zarr = isinstance(output, zarr.hierarchy.Group)
    for name in adata.obsm.keys():
        m = adata.obsm[name]
        dimensions = 2
        if m.shape[1] > 2:
            dimensions = 3
        d = {}
        for i in range(dimensions):
            d[name + '_' + str(i + 1)] = m[:, i].astype('float32')
        if is_zarr:
            g = output.create_group('/{}'.format(name))
            for key in d:
                value = d[key]
                ds = g.create_dataset(key, shape=value.shape, chunks=None, dtype=value.dtype)
                ds[:] = value
        else:
            write_table(d, output, name)


def save_data_obs(adata, output):
    logger.info('writing adata obs')
    is_zarr = isinstance(output, zarr.hierarchy.Group)
    for name in adata.obs:
        # TODO sort?, row group size?
        value = adata.obs[name]
        if is_zarr:
            g = output.create_group('/{}'.format(name))
            object_codec = None
            dtype = value.dtype
            if pd.api.types.is_categorical_dtype(value):
                object_codec = numcodecs.Categorize(value.cat.categories, dtype=object)
                dtype = object
                value = value.astype(str)

            value = value.values
            ds = g.create_dataset('value', shape=value.shape, chunks=None, dtype=dtype, object_codec=object_codec)
            ds[:] = value
        else:
            write_table(dict(value=value), output, name)
    value = adata.obs.index.values
    if is_zarr:
        g = output.create_group('/{}'.format('index'))
        ds = g.create_dataset('value', shape=value.shape, chunks=None, dtype=value.dtype,
            object_codec=numcodecs.MsgPack())
        ds[:] = value
    else:
        write_table(dict(value=value), output, 'index')


def save_adata_X_chunk(adata, adata_col_slice, output):
    X_slice = adata.X[:, adata_col_slice]
    names = adata.var.index[adata_col_slice]
    is_sparse = scipy.sparse.issparse(X_slice)
    if is_sparse and scipy.sparse.isspmatrix_csr(X_slice):
        X_slice = X_slice.tocsc()
    is_zarr = isinstance(output, zarr.hierarchy.Group)
    for j in range(X_slice.shape[1]):
        X = X_slice[:, j]
        if is_sparse:
            X = X.toarray().flatten()
        indices = np.where(X != 0)[0]
        values = X[indices]
        if is_zarr:
            g = output.create_group('/{}'.format(names[j]))
            ds = g.create_dataset('index', shape=indices.shape, chunks=None, dtype=indices.dtype)
            ds[:] = indices
            ds = g.create_dataset('value', shape=values.shape, chunks=None, dtype=values.dtype)
            ds[:] = values
        else:
            write_table(dict(index=indices, value=values), output, names[j])
        if j > 0 and (j + 1) % 1000 == 0:
            logger.info('Wrote adata X {}/{}'.format(j + 1, X_slice.shape[1]))
