import numpy as np
import pandas as pd
import pytest

from cirrocumulus.data_processing import handle_embedding
from cirrocumulus.embedding_aggregator import EmbeddingAggregator, get_basis


def create_df(test_data, measures, dimensions, basis):
    df = pd.DataFrame(test_data[:, measures].X.toarray(), columns=measures)
    df = df.join(test_data.obs[dimensions].reset_index())
    return df.join(pd.DataFrame(test_data.obsm['X_umap'][:, 0:2], columns=basis['coordinate_columns']))


def test_no_binning(dataset_api, input_dataset, test_data, measures, dimensions, continuous_obs, basis):
    basis = get_basis(basis)
    results = handle_embedding(dataset_api=dataset_api, dataset=input_dataset, basis=basis,
        measures=measures + list(map(lambda x: 'obs/' + x, continuous_obs)), dimensions=dimensions)

    df = create_df(test_data, measures, dimensions + continuous_obs, basis)
    for key in basis['coordinate_columns']:
        np.testing.assert_array_equal(results['coordinates'][key], df[key].values, err_msg=key)
    for key in measures:
        np.testing.assert_array_equal(results['values'][key], df[key].values, err_msg=key)
    for key in dimensions:
        np.testing.assert_array_equal(results['values'][key], test_data.obs[key].values, err_msg=key)


def group_df(test_data, measures, dimensions, continuous_obs, basis):
    agg_dict = {}
    df = create_df(test_data, measures, dimensions + continuous_obs, basis)
    for key in basis['coordinate_columns']:
        agg_dict[key] = 'max'
    for key in measures + continuous_obs:
        agg_dict[key] = 'max'
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
        np.testing.assert_array_equal(results['values'][key], grouped_df[key].values, err_msg=key)
    for key in continuous_obs:
        np.testing.assert_array_equal(results['values'][key], grouped_df[key].values, err_msg=key)
    for key in dimensions:
        np.testing.assert_array_equal(results['values'][key]['value'], grouped_df[key].values, err_msg=key)


@pytest.fixture(autouse=True)
def quick():
    return [True, False]


def test_binning(dataset_api, input_dataset, test_data, measures, dimensions, continuous_obs, basis, quick):
    nbins = 100
    agg_function = 'max'
    basis = get_basis(basis, nbins, agg_function)
    grouped_df = group_df(test_data, measures, dimensions, continuous_obs, basis)
    results = handle_embedding(dataset_api=dataset_api, dataset=input_dataset,
        basis=basis, measures=measures + list(map(lambda x: 'obs/' + x, continuous_obs)), dimensions=dimensions,
        quick=quick)
    diff_binning(grouped_df, measures, dimensions, continuous_obs, basis, results)


# def test_saved_embedding(tmp_path, dataset_api, test_data, measures, dimensions, continuous_obs, basis):
#     import zarr
#     url = str(tmp_path) + '.zarr'
#     make_ordered(test_data.obs, None)
#     test_data = test_data[:, measures]
#     test_data.write_zarr(url)
#     store = zarr.open(url)
#     nbins = 100
#     agg_function = 'max'
#     basis_name = basis
#     obsm = test_data.obsm[basis_name]
#     basis = get_basis(basis, nbins, agg_function)
#     for i in range(len(basis['coordinate_columns'])):
#         test_data.obs[basis['coordinate_columns'][i]] = obsm[:, i]
#
#     EmbeddingAggregator.convert_coords_to_bin(df=test_data.obs,
#         nbins=basis['nbins'],
#         coordinate_columns=basis['coordinate_columns'],
#         bin_name=basis['full_name'])
#
#     result = EmbeddingAggregator(obs_measures=continuous_obs,
#         var_measures=[], dimensions=dimensions,
#         count=False,
#         nbins=basis['nbins'], basis=basis, agg_function=basis['agg']).execute(test_data)
#     basis_group = require_binned_basis_group(store, basis)
#     obs_group = basis_group.require_group('obs')
#     coords_group = basis_group.require_group('coords')
#     write_basis_obs(basis, coords_group, obs_group, result)
#     result = EmbeddingAggregator(obs_measures=[],
#         var_measures=measures, dimensions=[],
#         count=False,
#         nbins=basis['nbins'], basis=basis, agg_function=basis['agg']).execute(test_data)
#     write_basis_X(coords_group, basis_group, test_data, result)
#
#     input_dataset = Entity(url, {'name': tmp_path, 'url': url})
#     result = dataset_api.read_precomputed_basis(input_dataset, obs_keys=dimensions + continuous_obs, var_keys=measures,
#         basis=basis)
#     grouped_df = group_df(test_data, measures, dimensions, continuous_obs, basis)
#     diff_binning(grouped_df, measures, dimensions, continuous_obs, basis, result)
