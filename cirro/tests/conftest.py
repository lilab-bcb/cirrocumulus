import os

import anndata
import pytest

from cirro.dataset_api import DatasetAPI
from cirro.entity import Entity
from cirro.h5ad_dataset import H5ADDataset


@pytest.fixture(scope='module', autouse=True)
def test_data():
    return anndata.read('test-data/pbmc3k_no_raw.h5ad')


@pytest.fixture(scope='module', autouse=True)
def measures():
    return ['TNFRSF4', 'DSCR3', 'SUMO3']


@pytest.fixture(scope='module', autouse=True)
def by():
    return 'louvain'


@pytest.fixture(scope='module', autouse=True)
def dimensions():
    return ['louvain']


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
    dataset_api.add(H5ADDataset(backed=h5_dataset_backed, force_sparse=h5_dataset_force_sparse))
    return dataset_api


@pytest.fixture(scope='module', params=['test-data/pbmc3k_no_raw.h5ad'])
def input_dataset(request):
    meta = {'name': os.path.splitext(os.path.basename(request.param))[0], 'url': request.param}
    if request.param.endswith('.json'):
        import json
        with open(request.param, 'rt') as f:
            meta.update(json.load(f))
    return Entity(request.param, meta)
