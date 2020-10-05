import logging
import os

import numpy as np
import pandas._libs.json as ujson
import scipy.sparse

logger = logging.getLogger("cirro")


def write_json(d, output_dir, name):
    os.makedirs(output_dir, exist_ok=True)
    with open(output_dir + os.path.sep + name + '.json', 'wt') as f:
        c = ujson.dumps(d, double_precision=2, orient='values')
        f.write(c)


def save_adata_json(adata, schema, output_directory):
    logger.info('Save adata')
    os.makedirs(output_directory, exist_ok=True)
    with open(os.path.join(output_directory, 'schema.json'), 'wt') as f:
        # json.dump(result, f)
        f.write(ujson.dumps(schema, double_precision=2, orient='values'))

    save_adata_X(adata, output_directory)
    save_data_obs(adata, output_directory)
    save_data_obsm(adata, output_directory)


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
        write_json(dict(index=indices, value=values), X_dir, names[j])
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
        write_json(d, obsm_dir, name)


def save_data_obs(adata, obs_dir):
    logger.info('writing adata obs')
    for name in adata.obs:
        # TODO sort?
        value = adata.obs[name]
        write_json(value, obs_dir, name)
    write_json(adata.obs.index.values, obs_dir, 'index')
