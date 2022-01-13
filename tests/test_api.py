import os

import anndata
import pytest

from cirrocumulus.envir import CIRRO_DB_URI, CIRRO_TEST
from cirrocumulus.launch import create_app, configure_app
from cirrocumulus.prepare_data import PrepareData
from cirrocumulus.serve import cached_app


@pytest.fixture(scope='session', params=[True, False])
def app_conf(request, tmpdir_factory):
    dataset_path = 'test-data/pbmc3k_no_raw.h5ad'

    if not request.param:
        app = create_app()
        configure_app(app, [dataset_path], None, None)
        dataset_id = dataset_path
        os.environ[CIRRO_TEST] = 'false'
    else:
        os.environ[CIRRO_TEST] = 'true'
        os.environ[CIRRO_DB_URI] = 'mongodb://localhost:27018/cirrocumulus-test'
        app = cached_app()
    with app.test_client() as client:
        if request.param:
            # insert dataset
            output_dir = str(tmpdir_factory.mktemp("data").join("test.zarr"))
            PrepareData(datasets=[anndata.read(dataset_path)], output=output_dir, output_format='zarr',
                        no_auto_groups=True).execute()
            r = client.post('/api/dataset', data=dict(url=output_dir, name='test'))
            dataset_id = r.json['id']
        yield client, dataset_id


def test_server(app_conf):
    client, dataset_id = app_conf
    r = client.get('/api/server').json
    assert r['email'] is None


def test_schema(app_conf):
    client, dataset_id = app_conf
    r = client.get('/api/schema?id={}'.format(dataset_id)).json
    assert isinstance(r['var'], list)
    assert isinstance(r['obs'], list)
    assert isinstance(r['embeddings'], list)
    assert len(r['obsCat']) == 1 and r['obsCat'][0] == 'louvain'
    assert r['shape'][0] == 2638 and r['shape'][1] == 1838
