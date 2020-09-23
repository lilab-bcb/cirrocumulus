import gzip
import json
import logging
import os

import numpy as np
import pandas as pd
import pandas._libs.json as ujson
import scipy.sparse

from cirrocumulus.simple_data import SimpleData

logger = logging.getLogger("cirro")


def write_json(d, f, name, index, compress):
    c = ujson.dumps(d, double_precision=2, orient='values').encode('UTF-8')
    if compress:
        c = gzip.compress(c)
    start = f.tell()
    end = start + len(c)
    index[name] = [start, end - 1]
    f.write(c)


def read_adata_json(path, keys):
    with gzip.open(path + '.idx.json', 'rt') as f:
        index = json.load(f)

    def read_key(key):
        start, end = index[key]
        f.seek(start)
        b = f.read(end - start)
        if xx:
            b = gzip.decompress(b)
        return json.loads(b.decode('UTF-8'))

    df = pd.DataFrame()
    with open(path, 'rb') as f:
        index = schema['index']
        shape = schema['shape']
        for key in keys:
            value = read_key(key)
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


def save_adata_json(adata, output_path, compress):
    logger.info('Save adata')

    index = {}  # key to byte start-end
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    data_file = os.path.join(output_path, 'data.json')
    with open(data_file, 'wb') as f:
        save_adata_X(adata, f, index, compress)
        save_data_obs(adata, f, index, compress)
        save_data_obsm(adata, f, index, compress)

    with open(os.path.join(output_path, 'index.json'), 'wt') as f:
        # json.dump(result, f)
        schema = SimpleData.schema(adata)
        schema['format'] = 'BGZF' if compress else 'JSON'
        result = dict(index=index, schema=schema, file=os.path.basename(data_file))
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
        write_json(dict(index=indices, value=values), f, names[j], index, compress)
        if j > 0 and (j + 1) % 1000 == 0:
            logger.info('Wrote adata X {}/{}'.format(j + 1, adata_X.shape[1]))


def save_data_obsm(adata, f, index, compress):
    logger.info('writing adata obsm')

    for name in adata.obsm.keys():
        m = adata.obsm[name]
        dimensions = 3 if m.shape[1] > 2 else 2
        d = {}
        for i in range(dimensions):
            d[name + '_' + str(i + 1)] = m[:, i].astype('float32')
        write_json(d, f, name, index, compress)


def save_data_obs(adata, f, index, compress):
    logger.info('writing adata obs')
    for name in adata.obs:
        value = adata.obs[name]
        write_json(value, f, name, index, compress)
    write_json(adata.obs.index.values, f, 'index', index, compress)
