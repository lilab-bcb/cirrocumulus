import anndata
import fsspec
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
import scipy.sparse
from cirrocumulus.anndata_util import get_base
from cirrocumulus.diff_exp import DE
from cirrocumulus.parquet_dataset import ParquetDataset
from cirrocumulus.prepare_data import PrepareData
from cirrocumulus.zarr_dataset import ZarrDataset
from numpy.random import negative_binomial, binomial, seed
from scipy import sparse as sp


def get_example_data(sparse=False):
    seed(1234)
    # create test object
    adata = anndata.AnnData(
        np.multiply(binomial(1, 0.15, (100, 20)), negative_binomial(2, 0.25, (100, 20)))
    )
    # adapt marker_genes for cluster (so as to have some form of reasonable input
    adata.X[0:10, 0:5] = np.multiply(
        binomial(1, 0.9, (10, 5)), negative_binomial(1, 0.5, (10, 5))
    )

    # The following construction is inefficient, but makes sure that the same data is used in the sparse case
    if sparse:
        adata.X = sp.csr_matrix(adata.X)

    # Create cluster according to groups
    adata.obs['sc_groups'] = pd.Categorical(
        np.concatenate(
            (
                np.zeros((10,), dtype=int),
                np.ones((90,), dtype=int),
            )
        )
    )

    return adata


@pytest.fixture(autouse=True, params=[True, False])
def sparse(request):
    return request.param


def diff_results(adata, obs_field, results, group='0'):
    sc.tl.rank_genes_groups(adata, obs_field, method='t-test', pts=True)
    rank_genes_groups = adata.uns['rank_genes_groups']
    sc_scores = rank_genes_groups['scores'][group]
    sc_pvals = rank_genes_groups['pvals'][group]
    sc_pts = rank_genes_groups['pts'][group]
    sc_lfc = rank_genes_groups['logfoldchanges'][group]

    sc_df = pd.DataFrame(index=rank_genes_groups['names'][group],
                         data={'pvals': sc_pvals, 'scores': sc_scores, 'pts': sc_pts, 'lfc': sc_lfc})
    sc_df = sc_df.loc[adata.var.index]
    np.testing.assert_allclose(sc_df['pvals'], results['pvals'])

    np.testing.assert_allclose(sc_df['scores'], results['scores'], atol=1e-015)
    np.testing.assert_allclose(sc_df['lfc'], results['logfoldchanges'], atol=1e-015)

    if scipy.sparse.issparse(adata.X):
        np.testing.assert_allclose(sc_df['pts'], results['frac_expressed1'])
    else:
        assert results['frac_expressed1'] is None


@pytest.mark.parametrize("file_format", ['zarr', 'parquet'])
def test_de_backed(sparse, file_format, tmp_path):
    fs = fsspec.filesystem('file')
    adata = get_example_data(sparse)
    output_dir = str(tmp_path)
    prepare_data = PrepareData(datasets=[adata], output=output_dir, output_format=file_format)
    prepare_data.execute()
    if file_format == 'parquet':
        reader = ParquetDataset()
    elif file_format == 'zarr':
        reader = ZarrDataset()
    batch_size = 30
    obs_field = 'sc_groups'
    nfeatures = adata.shape[1]

    def get_batch_fn(i):
        end = min(nfeatures, i + batch_size)
        return reader.read_dataset(filesystem=fs, path=output_dir, dataset=dict(id=''),
                                   keys=dict(X=[slice(i, end)]))

    results = DE(series=adata.obs[obs_field], nfeatures=nfeatures, batch_size=batch_size, get_batch_fn=get_batch_fn,
                 base=get_base(adata), one_vs_rest=True)
    diff_results(adata, obs_field, results.pair2results[0])


def test_de_2_groups(sparse):
    adata = get_example_data(sparse)
    batch_size = 3
    obs_field = 'sc_groups'
    nfeatures = adata.shape[1]
    get_batch_fn = lambda i: adata[:, i:min(nfeatures, i + batch_size)]

    results = DE(series=adata.obs[obs_field], nfeatures=nfeatures, batch_size=batch_size, get_batch_fn=get_batch_fn,
                 base=get_base(adata), one_vs_rest=True)
    diff_results(adata, obs_field, results.pair2results[0])


def test_de_4_groups(sparse):
    adata1 = get_example_data(sparse)
    adata2 = get_example_data(sparse)
    adata2.obs['sc_groups'] = adata2.obs['sc_groups'].replace({0: 2, 1: 3})
    adata = anndata.concat((adata1, adata2))
    adata.obs_names_make_unique()
    batch_size = 3
    obs_field = 'sc_groups'
    adata.obs[obs_field] = adata.obs[obs_field].astype('category')
    nfeatures = adata.shape[1]
    get_batch_fn = lambda i: adata[:, i:min(nfeatures, i + batch_size)]
    de = DE(series=adata.obs[obs_field], nfeatures=nfeatures, batch_size=batch_size, get_batch_fn=get_batch_fn,
            base=get_base(adata))
    for i in range(4):
        diff_results(adata, obs_field, de.pair2results[i], str(i))
