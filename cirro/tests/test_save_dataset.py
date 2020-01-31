import shutil

import anndata
import numpy as np
import scipy.sparse

from cirro.dataset_api import DatasetAPI
from cirro.embedding_aggregator import get_basis
from cirro.entity import Entity
from cirro.io import save_adata
from cirro.parquet_dataset import ParquetDataset


def test_save_adata(measures, dimensions, basis, h5_dataset_force_sparse):
    dataset_api = DatasetAPI()
    dataset_api.add(ParquetDataset())
    test_data = anndata.read('test-data/pbmc3k_no_raw.h5ad')

    test_data = test_data[:, measures].copy()
    if h5_dataset_force_sparse and not scipy.sparse.issparse(test_data.X):
        test_data.X = scipy.sparse.csr_matrix(test_data.X)
    save_adata(test_data, 'test_adata/data')

    # with open('test_adata/index.json', 'wt') as f:
    #     json.dump(SimpleData.schema(test_data), f)
    basis = get_basis(basis, -1, None)
    meta = {'name': 'test_adata/data', 'url': 'test_adata/index.json', 'nObs': test_data.shape[0]}
    input_dataset = Entity('test_adata/index.json', meta)
    result = dataset_api.read(input_dataset, obs_keys=dimensions, var_keys=measures, basis=[basis])
    result.X = result.X.toarray()
    if scipy.sparse.issparse(test_data.X):
        test_data.X = test_data.X.toarray()
    shutil.rmtree('test_adata')
    for dimension in dimensions:
        assert (test_data.obs[dimension].values != result.obs[dimension].values).sum() == 0
    np.testing.assert_array_almost_equal(result.X, test_data.X)
