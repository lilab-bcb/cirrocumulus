import numpy as np
import pandas as pd

from cirro.data_processing import process_data
from cirro.embedding_aggregator import EmbeddingAggregator


def create_df(test_data, measures, dimensions, basis):
    df = pd.DataFrame(test_data[:, measures].X.toarray(), columns=measures)
    df = df.join(test_data.obs[dimensions].reset_index())
    return df.join(pd.DataFrame(test_data.obsm['X_umap'][:, 0:2], columns=basis['coordinate_columns']))


def test_no_binning(dataset_api, input_dataset, test_data, measures, dimensions, basis):
    process_results = process_data(dataset_api=dataset_api, dataset=input_dataset, basis=basis,
        embedding_measures=measures, embedding_dimensions=dimensions,
        return_types=['embedding'])

    results = process_results['embedding']
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
        coordinate_column_to_range=None)
    # group by bin then agg

    grouped_df = df.groupby(df.index).agg(agg_dict)

    process_results = process_data(dataset_api=dataset_api, dataset=input_dataset, basis=basis, nbins=nbins,
        embedding_measures=measures, embedding_dimensions=dimensions,
        return_types=['embedding'], agg_function=agg_function)

    results = process_results['embedding']
    for key in basis['coordinate_columns']:
        np.testing.assert_array_equal(results['coordinates'][key], grouped_df[key].values, err_msg=key)
    for key in measures:
        np.testing.assert_array_equal(results['values'][key], grouped_df[key].values, err_msg=key)
    for key in dimensions:
        np.testing.assert_array_equal(results['values'][key]['mode'], grouped_df[key].values, err_msg=key)
