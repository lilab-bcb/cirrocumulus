import logging
import os

import anndata
import numpy as np
import pytest
import scipy.sparse

from cirro.dataset_api import DatasetAPI
from cirro.embedding_aggregator import get_basis
from cirro.entity import Entity
from cirro.io import save_adata
from cirro.parquet_dataset import ParquetDataset
from cirro.zarr_dataset import ZarrDataset

logger = logging.getLogger("cirro")


@pytest.fixture(params=['parquet', 'zarr'])
def output_format(request):
    return request.param


def test_save_adata(tmp_path, measures, dimensions, basis, h5_dataset_force_sparse, output_format):
    dataset_api = DatasetAPI()
    dataset_api.add(ParquetDataset())
    dataset_api.add(ZarrDataset())
    data_path = os.path.join(str(tmp_path), 'data')
    test_data = anndata.read('test-data/pbmc3k_no_raw.h5ad')
    test_data = test_data[:, measures].copy()
    if h5_dataset_force_sparse and not scipy.sparse.issparse(test_data.X):
        test_data.X = scipy.sparse.csr_matrix(test_data.X)
    save_adata(test_data, data_path, output_format=output_format)

    extension = 'pjson' if output_format == 'parquet' else 'zjson'
    basis = get_basis(basis)

    url = os.path.join(tmp_path, 'index.{}'.format(extension))
    meta = {'name': tmp_path, 'url': url, 'nObs': test_data.shape[0]}
    input_dataset = Entity(url, meta)
    result = dataset_api.read(input_dataset, obs_keys=dimensions, var_keys=measures, basis=[basis])
    result.X = result.X.toarray()
    if scipy.sparse.issparse(test_data.X):
        test_data.X = test_data.X.toarray()
    for dimension in dimensions:
        assert (test_data.obs[dimension].values != result.obs[dimension].values).sum() == 0
    np.testing.assert_array_almost_equal(result.X, test_data.X)
