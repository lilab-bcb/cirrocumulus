import numpy as np
import pandas as pd
import pytest

from cirrocumulus.data_processing import handle_data
from cirrocumulus.embedding_aggregator import EmbeddingAggregator, get_basis


def create_df(test_data, measures, dimensions, basis):
    df = pd.DataFrame(test_data[:, measures].X.toarray(), columns=measures)
    df = df.join(test_data.obs[dimensions].reset_index())
    return df.join(pd.DataFrame(test_data.obsm['X_umap'][:, 0:2], columns=basis['coordinate_columns']))


def group_df(test_data, measures, dimensions, continuous_obs, basis):
    df = create_df(test_data, measures, dimensions + continuous_obs, basis)
    agg_dict = {}
    for key in basis['coordinate_columns']:
        agg_dict[key] = 'max'
    for key in measures + continuous_obs:
        agg_dict[key] = basis['agg']
    for key in dimensions:
        agg_dict[key] = lambda x: x.mode()[0]
    EmbeddingAggregator.convert_coords_to_bin(df, nbins=100,
        coordinate_columns=basis['coordinate_columns'],
        bin_name=basis['full_name'])
    # group by bin then agg
    return df.groupby(basis['full_name']).agg(agg_dict)


def diff_binning(grouped_df, measures, dimensions, continuous_obs, basis, results):
    for key in basis['coordinate_columns']:
        np.testing.assert_array_equal(results['coordinates'][key], grouped_df[key].values, err_msg=key)
    for key in measures:
        np.testing.assert_allclose(results['values'][key], grouped_df[key].values, err_msg=key, rtol=1e-05, atol=1e-06)
    for key in continuous_obs:
        np.testing.assert_allclose(results['values'][key], grouped_df[key].values, err_msg=key, rtol=1e-05,
            atol=1e-08)
    for key in dimensions:
        np.testing.assert_array_equal(results['values'][key]['value'], grouped_df[key].values, err_msg=key)


def test_no_binning(dataset_api, input_dataset, test_data, measures, dimensions, continuous_obs, basis):
    embedding_list = [dict(ndim=2, basis=basis)]
    basis = get_basis(basis)
    values = dict(dimensions=dimensions, measures=measures + list(map(lambda x: 'obs/' + x, continuous_obs)))
    results = handle_data(dataset_api=dataset_api, dataset=input_dataset, values=values, embedding_list=embedding_list)
    df = create_df(test_data, measures, dimensions + continuous_obs, basis)
    for key in basis['coordinate_columns']:
        np.testing.assert_array_equal(results['embeddings'][0]['coordinates'][key], df[key].values, err_msg=key)
    for key in measures:
        np.testing.assert_array_equal(results['values'][key], df[key].values, err_msg=key)
    for key in dimensions:
        np.testing.assert_array_equal(results['values'][key], test_data.obs[key], err_msg=key)


@pytest.fixture(autouse=True, params=['sum', 'mean', 'max'])
def agg_function(request):
    return request.param


def test_binning(dataset_api, input_dataset, test_data, measures, dimensions, continuous_obs, basis,
                 agg_function):
    nbins = 100
    basis_dict = get_basis(basis, nbins, agg_function)
    grouped_df = group_df(test_data, measures, dimensions, continuous_obs, basis_dict)
    embedding_list = [dict(nbins=nbins, ndim=2, basis=basis, dimensions=dimensions, agg=agg_function)]
    embedding_list[0]['measures'] = measures + list(map(lambda x: 'obs/' + x, continuous_obs))
    results = handle_data(dataset_api=dataset_api, dataset=input_dataset, embedding_list=embedding_list)
    diff_binning(grouped_df, measures, dimensions, continuous_obs, basis_dict, results['embeddings'][0])
