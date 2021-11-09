import pytest

from cirrocumulus.launch import create_app, configure_app


@pytest.fixture
def client():
    app = create_app()
    configure_app(app, ['test-data/pbmc3k_no_raw.h5ad'], None, None)
    with app.test_client() as client:
        yield client


def test_server(client):
    r = client.get('/api/server').json
    assert r['email'] is None


def test_schema(client):
    r = client.get('/api/schema?id=test-data/pbmc3k_no_raw.h5ad').json
    assert isinstance(r['var'], list)
    assert isinstance(r['obs'], list)
    assert isinstance(r['embeddings'], list)
    assert len(r['obsCat']) == 1 and r['obsCat'][0] == 'louvain'
    assert r['shape'][0] == 2638 and r['shape'][1] == 1838
