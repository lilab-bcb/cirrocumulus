import os

import anndata
import pytest
import scipy.sparse

from cirrocumulus.anndata_dataset import AnndataDataset
from cirrocumulus.dataset_api import DatasetAPI


@pytest.fixture(scope='module', autouse=True, params=[True, False])
def test_data(request):
    adata = anndata.read('test-data/pbmc3k_no_raw.h5ad')
    if request.param:
        adata.X = scipy.sparse.csr_matrix(adata.X)
    return adata


@pytest.fixture(scope='module', autouse=True, params=['small', 'large'])
def measures(request):
    return ['DSCR3', 'TNFRSF4', 'SUMO3'] if request.param == 'small' else list(
        anndata.read('test-data/pbmc3k_no_raw.h5ad').var.index[0:50])


@pytest.fixture(scope='module', autouse=True)
def by():
    return 'louvain'


@pytest.fixture(scope='module', autouse=True)
def dimensions():
    return ['louvain']


@pytest.fixture(scope='module', autouse=True)
def continuous_obs():
    return ['n_genes', 'percent_mito']


@pytest.fixture(scope='module')
def basis():
    return 'X_umap'


@pytest.fixture(scope='module', params=['r', None])
def h5_dataset_backed(request):
    return request.param


@pytest.fixture(scope='module', params=[True, False])
def h5_dataset_force_sparse(request):
    return request.param


@pytest.fixture(scope='module', autouse=True)
def dataset_api(h5_dataset_backed):
    dataset_api = DatasetAPI()
    dataset_api.add(AnndataDataset(backed=h5_dataset_backed))
    return dataset_api


@pytest.fixture(scope='module', params=['test-data/pbmc3k_no_raw.h5ad'])
def input_dataset(request):
    meta = {'name': os.path.splitext(os.path.basename(request.param))[0], 'url': request.param, 'id': request.param}
    if request.param.endswith('.json'):
        import json
        with open(request.param, 'rt') as f:
            meta.update(json.load(f))

    return meta
