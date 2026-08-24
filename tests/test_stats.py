import numpy as np
import pytest
import scipy.sparse

from cirrocumulus.anndata_util import X_stats
from cirrocumulus.feature_aggregator import FeatureAggregator


@pytest.fixture(params=["sparse", "dense"])
def stats_data(request, test_data):
    adata = test_data.copy()
    if request.param == "sparse":
        adata.X = scipy.sparse.csr_matrix(adata.X)
    elif scipy.sparse.issparse(adata.X):
        adata.X = adata.X.toarray()
    return adata


def test_X_stats_is_one_dimensional(stats_data):
    df = X_stats(stats_data)
    assert list(df.index) == list(stats_data.var.index)
    for column in ("min", "max", "sum", "mean", "numExpressed"):
        assert df[column].values.ndim == 1


def test_X_stats_matches_dense(stats_data):
    df = X_stats(stats_data)
    X = stats_data.X
    dense = X.toarray() if scipy.sparse.issparse(X) else np.asarray(X)
    # sparse and dense reductions accumulate in a different order, so float32 sums of
    # near-zero scaled values differ in the last bits -- compare with an absolute floor
    scale = float(np.abs(dense).max())
    np.testing.assert_allclose(df["min"].values, dense.min(axis=0), rtol=1e-5, atol=1e-5 * scale)
    np.testing.assert_allclose(df["max"].values, dense.max(axis=0), rtol=1e-5, atol=1e-5 * scale)
    np.testing.assert_allclose(
        df["sum"].values, dense.sum(axis=0), rtol=1e-4, atol=1e-4 * scale * dense.shape[0] ** 0.5
    )
    np.testing.assert_allclose(df["mean"].values, dense.mean(axis=0), rtol=1e-4, atol=1e-5 * scale)
    np.testing.assert_array_equal(df["numExpressed"].values, np.count_nonzero(dense, axis=0))


def test_feature_aggregator_summary(stats_data, measures, dimensions):
    result = FeatureAggregator(var_measures=measures, dimensions=dimensions).execute(
        stats_data[:, measures]
    )
    for dimension in dimensions:
        counts = result[dimension]["counts"]
        assert sum(counts) == stats_data.shape[0]
    for measure in measures:
        assert set(result[measure]) == {"min", "max", "sum", "mean", "numExpressed"}
