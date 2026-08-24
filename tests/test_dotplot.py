import numpy as np
import pandas as pd
import pytest
import scipy.sparse

from cirrocumulus.dotplot_aggregator import DotPlotAggregator


@pytest.fixture(params=["sparse", "dense"])
def dotplot_data(request, test_data):
    adata = test_data.copy()
    if request.param == "sparse":
        adata.X = scipy.sparse.csr_matrix(adata.X)
    elif scipy.sparse.issparse(adata.X):
        adata.X = adata.X.toarray()
    return adata


def _expected(adata, measure, by):
    values = adata[:, measure].X
    values = np.asarray(values.todense() if scipy.sparse.issparse(values) else values).ravel()
    grouped = pd.Series(values, index=adata.obs.index).groupby(adata.obs[by].values, observed=True)
    return grouped.mean(), grouped.apply(lambda x: 100 * (x.values != 0).sum() / len(x))


def test_dotplot_accepts_anndata(dotplot_data, measures, by):
    results = DotPlotAggregator(var_measures=measures, dimensions=[by]).execute(dotplot_data)
    assert len(results) == 1
    result = results[0]
    assert result["name"] == by
    assert list(result["categories"]) == list(dotplot_data.obs[by].cat.categories)
    assert [v["name"] for v in result["values"]] == list(measures)
    for value in result["values"]:
        assert len(value["mean"]) == len(result["categories"])
        assert np.all(np.asarray(value["percentExpressed"]) >= 0)
        assert np.all(np.asarray(value["percentExpressed"]) <= 100)


def test_dotplot_values_match_manual_groupby(dotplot_data, measures, by):
    result = DotPlotAggregator(var_measures=measures, dimensions=[by]).execute(dotplot_data)[0]
    categories = list(result["categories"])
    for value in result["values"]:
        expected_mean, expected_pct = _expected(dotplot_data, value["name"], by)
        np.testing.assert_allclose(
            np.asarray(value["mean"]), expected_mean.loc[categories].values, rtol=1e-4, atol=1e-5
        )
        np.testing.assert_allclose(
            np.asarray(value["percentExpressed"]), expected_pct.loc[categories].values, rtol=1e-4
        )


def test_dotplot_multiple_dimensions(dotplot_data, measures, by):
    adata = dotplot_data
    adata.obs["half"] = pd.Categorical(np.where(np.arange(adata.shape[0]) % 2 == 0, "a", "b"))
    results = DotPlotAggregator(var_measures=measures, dimensions=[[by, "half"]]).execute(adata)
    assert len(results) == 1
    assert results[0]["name"] == "{}-half".format(by)
