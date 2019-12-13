import numpy as np
import pandas as pd

from cirro.data_processing import process_data
from cirro.embedding_aggregator import EmbeddingAggregator


def test_dimension_filter(dataset_api, input_dataset, test_data):
    process_results = process_data(dataset_api=dataset_api, dataset=input_dataset, return_types=['ids'],
        data_filter={'filters': [['obs/louvain', 'in', ['1', '5']]]})
    ids_summary = process_results['ids']
    matched_data = test_data[test_data.obs['louvain'].isin(['5', '1'])]
    assert len(np.intersect1d(ids_summary, matched_data.obs.index)) == matched_data.shape[0], '{}, {}'.format(
        len(ids_summary), matched_data.shape[0])


def test_binning_dimension_filter(dataset_api, input_dataset, test_data, measures, dimensions, basis):
    nbins = 100
    df = pd.DataFrame(test_data[:, measures].X.toarray(), columns=measures)
    df = df.join(test_data.obs[dimensions])
    df['id'] = test_data.obs.index
    df = df.join(pd.DataFrame(test_data.obsm['X_umap'][:, 0:2], columns=basis['coordinate_columns']))
    EmbeddingAggregator.convert_coords_to_bin(df, nbins=nbins,
        coordinate_columns=basis['coordinate_columns'],
        coordinate_column_to_range=None)
    df = df[df['louvain'].isin(['5', '1'])]
    process_results = process_data(dataset_api=dataset_api, basis=basis, nbins=nbins, dataset=input_dataset,
        return_types=['ids'], data_filter={'filters': [['obs/louvain', 'in', ['1', '5']]]})
    ids_summary = process_results['ids']
    both = np.intersect1d(ids_summary, df['id'])
    assert len(both) == df.shape[0], '{}, {}'.format(len(ids_summary), len(df))


def test_measure_filter(dataset_api, input_dataset, test_data):
    process_results = process_data(dataset_api=dataset_api, dataset=input_dataset, return_types=['ids'],
        data_filter={'filters': [['DSCR3', '>', 2]]})
    ids_summary = process_results['ids']
    matched_data = test_data[test_data[:, 'DSCR3'].X > 2]
    assert len(np.intersect1d(ids_summary, matched_data.obs.index)) == matched_data.shape[0]


def test_binning_measure_filter(dataset_api, input_dataset, test_data, measures, dimensions, basis):
    nbins = 100
    df = pd.DataFrame(test_data[:, measures].X.toarray(), columns=measures)
    df = df.join(test_data.obs[dimensions])
    df['id'] = test_data.obs.index
    df = df.join(pd.DataFrame(test_data.obsm['X_umap'][:, 0:2], columns=basis['coordinate_columns']))
    EmbeddingAggregator.convert_coords_to_bin(df, nbins=nbins,
        coordinate_columns=basis['coordinate_columns'],
        coordinate_column_to_range=None)
    df = df[df['DSCR3'] > 2]
    process_results = process_data(dataset_api=dataset_api, basis=basis, nbins=nbins, dataset=input_dataset,
        return_types=['ids'], data_filter={'filters': [['DSCR3', '>', 2]]})
    ids_summary = process_results['ids']
    both = np.intersect1d(ids_summary, df['id'])
    assert len(both), df.shape[0]
