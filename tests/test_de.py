import anndata
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from numpy.random import negative_binomial, binomial, seed
from scipy import sparse as sp

from cirrocumulus.de import DE


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


def test_de(sparse):
    adata = get_example_data(sparse)

    batch_size = 3
    obs_field = 'sc_groups'
    nfeatures = adata.shape[1]
    get_batch_fn = lambda i: adata[:, i:min(nfeatures, i + batch_size)]
    de = DE(adata, obs_field, nfeatures, batch_size, get_batch_fn)

    sc.tl.rank_genes_groups(adata, obs_field, method='t-test', pts=True)
    rank_genes_groups = adata.uns['rank_genes_groups']
    sc_scores = rank_genes_groups['scores']['0']
    sc_pvals = rank_genes_groups['pvals']['0']
    sc_pts = rank_genes_groups['pts']
    sc_lfc = rank_genes_groups['logfoldchanges']['0']

    features = rank_genes_groups['names']['0']
    sc_df = pd.DataFrame(index=features,
                         data={'pvals': sc_pvals, 'scores': sc_scores, 'pts': sc_pts['0'], 'pts_rest': sc_pts['1'],
                               'lfc': sc_lfc})
    sc_df = sc_df.loc[adata.var.index]
    np.testing.assert_allclose(sc_df['scores'], de.scores)
    np.testing.assert_allclose(sc_df['pvals'], de.pvals)
    np.testing.assert_allclose(sc_df['lfc'], de.logfoldchanges)
    if sparse:
        np.testing.assert_allclose(sc_df['pts'], de.frac_expressed.iloc[0])
        np.testing.assert_allclose(sc_df['pts_rest'], de.frac_expressed.iloc[1])
    else:
        assert de.frac_expressed is None

