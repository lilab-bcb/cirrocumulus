import numpy as np
import pandas as pd

from cirro.data_processing import handle_embedding
from cirro.embedding_aggregator import EmbeddingAggregator, get_basis


def create_df(test_data, measures, dimensions, basis):
    df = pd.DataFrame(test_data[:, measures].X.toarray(), columns=measures)
    df = df.join(test_data.obs[dimensions].reset_index())
    return df.join(pd.DataFrame(test_data.obsm['X_umap'][:, 0:2], columns=basis['coordinate_columns']))


def test_no_binning(dataset_api, input_dataset, test_data, measures, dimensions, basis):
    basis = get_basis(basis)
    results = handle_embedding(dataset_api=dataset_api, dataset=input_dataset, basis=basis,
        measures=measures, dimensions=dimensions)

    df = create_df(test_data, measures, dimensions, basis)
    for key in basis['coordinate_columns']:
        np.testing.assert_array_equal(results['coordinates'][key], df[key].values, err_msg=key)
    for key in measures:
        np.testing.assert_array_equal(results['values'][key], df[key].values, err_msg=key)
    for key in dimensions:
        np.testing.assert_array_equal(results['values'][key], test_data.obs[key].values, err_msg=key)


def test_binning(dataset_api, input_dataset, test_data, measures, dimensions, basis):
    nbins = 100
    agg_function = 'max'
    basis = get_basis(basis, nbins, agg_function)
    agg_dict = {}
    df = create_df(test_data, measures, dimensions, basis)
    for key in basis['coordinate_columns']:
        agg_dict[key] = 'max'
    for key in measures:
        agg_dict[key] = agg_function
    for key in dimensions:
        agg_dict[key] = lambda x: x.mode()[0]
    EmbeddingAggregator.convert_coords_to_bin(df, nbins=nbins,
        coordinate_columns=basis['coordinate_columns'],
        bin_name=basis['full_name'])
    # group by bin then agg

    grouped_df = df.groupby(basis['full_name']).agg(agg_dict)

    results = handle_embedding(dataset_api=dataset_api, dataset=input_dataset,
        basis=basis,
        measures=measures, dimensions=dimensions)

    for key in basis['coordinate_columns']:
        np.testing.assert_array_equal(results['coordinates'][key], grouped_df[key].values, err_msg=key)
    for key in measures:
        np.testing.assert_array_equal(results['values'][key], grouped_df[key].values, err_msg=key)
    for key in dimensions:
        np.testing.assert_array_equal(results['values'][key]['value'], grouped_df[key].values, err_msg=key)
