import os

import anndata
import pytest

from cirrocumulus.anndata_dataset import AnndataDataset
from cirrocumulus.dataset_api import DatasetAPI
from cirrocumulus.entity import Entity
from cirrocumulus.zarr_dataset_backed import ZarrDatasetBacked


@pytest.fixture(scope='module', autouse=True)
def test_data():
    return anndata.read('test-data/pbmc3k_no_raw.h5ad')


@pytest.fixture(scope='module', autouse=True, params=['small', 'large'])
def measures(request):
    return ['TNFRSF4', 'DSCR3', 'SUMO3'] if request.param == 'small' else list(
        anndata.read('test-data/pbmc3k_no_raw.h5ad').var.index)


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
def dataset_api(h5_dataset_backed, h5_dataset_force_sparse):
    dataset_api = DatasetAPI()
    dataset_api.add(AnndataDataset(backed=h5_dataset_backed, force_sparse=h5_dataset_force_sparse, extensions=['h5ad']))
    dataset_api.add(ZarrDatasetBacked())
    return dataset_api


@pytest.fixture(scope='module', params=['test-data/pbmc3k_no_raw.h5ad'])
def input_dataset(request):
    meta = {'name': os.path.splitext(os.path.basename(request.param))[0], 'url': request.param}
    if request.param.endswith('.json'):
        import json
        with open(request.param, 'rt') as f:
            meta.update(json.load(f))
    return Entity(request.param, meta)
