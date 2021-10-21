import gzip
import json
import logging
import os

import numpy as np
import pandas as pd
import pandas._libs.json as ujson
import scipy.sparse

logger = logging.getLogger("cirro")

LINE_END = '\n'.encode('UTF-8')


def write_jsonl(d, f, name, index, compress=False):
    output = {}
    output[name] = d
    c = ujson.dumps(output, double_precision=2, orient='values').encode('UTF-8')
    if compress:
        c = gzip.compress(c)
    start = f.tell()
    end = start + len(c)
    index[name] = [start, end - 1]
    f.write(c)
    f.write(LINE_END)


def read_adata_jsonl(path, keys):
    with open(path + '.idx.json', 'rt') as f:
        index = json.load(f)
    compress = False

    def read_key(f, index_key):
        start, end = index['index'][index_key]
        f.seek(start)
        b = f.read(end - start + 1)
        if compress:
            b = gzip.decompress(b)
        return json.loads(b.decode('UTF-8'))

    df = pd.DataFrame()
    with open(path, 'rb') as f:
        shape = read_key('schema')['shape']
        for key in keys:
            value = read_key(f, key)
            if isinstance(value, dict):
                if 'index' in value:
                    array = np.zeros(shape[0])
                    array[value['index']] = value['value']
                    df[key] = pd.arrays.SparseArray(array)
                else:
                    # obsm
                    for obsm_index in value:
                        df[obsm_index] = value[obsm_index]
            else:
                df[key] = value

    return df


def save_dataset_jsonl(dataset, schema, output_dir, base_name, filesystem):
    compress = False
    index = {}  # key to byte start-end
    filesystem.makedirs(output_dir, exist_ok=True)
    jsonl_path = os.path.join(output_dir, base_name)
    with filesystem.open(jsonl_path, 'wb') as f:
        save_adata_X(dataset, f, index, compress)
        save_data_obs(dataset, f, index, compress)
        save_data_obsm(dataset, f, index, compress)
        write_jsonl(schema, f, 'schema', index)

    with filesystem.open(os.path.join(output_dir, base_name + '.idx.json'), 'wt') as f:  # save index
        # json.dump(result, f)
        result = dict(index=index, file=os.path.basename(jsonl_path))
        f.write(ujson.dumps(result, double_precision=2, orient='values'))


def save_adata_X(adata, f, index, compress):
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
            write_jsonl(dict(index=indices, value=values), f, names[j], index, compress)
        else:
            write_jsonl(dict(value=X), f, names[j], index, compress)
        if j > 0 and (j + 1) % 1000 == 0:
            logger.info('Wrote adata X {}/{}'.format(j + 1, adata_X.shape[1]))


def save_data_obsm(adata, f, index, compress):
    logger.info('writing adata obsm')

    for name in adata.obsm.keys():
        value = adata.obsm[name]
        dim = value.shape[1]
        d = {}
        for i in range(dim):
            d[name + '_' + str(i + 1)] = value[:, i]
        write_jsonl(d, f, name, index, compress)


def save_data_obs(adata, f, index, compress):
    logger.info('writing adata obs')
    for name in adata.obs:
        series = adata.obs[name]
        value = series
        if pd.api.types.is_categorical_dtype(series):
            value = dict(values=series.values.codes, categories=series.cat.categories.values)
        write_jsonl(value, f, name, index, compress)
    write_jsonl(adata.obs.index.values, f, 'index', index, compress)
