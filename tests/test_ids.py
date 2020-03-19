import numpy as np

from cirrocumulus.data_processing import handle_selection_ids


def test_dimension_filter(dataset_api, input_dataset, test_data):
    process_results = handle_selection_ids(dataset_api=dataset_api, dataset=input_dataset,
        data_filter={'filters': [['obs/louvain', 'in', ['1', '5']]]})
    ids_summary = process_results['ids']
    matched_data = test_data[test_data.obs['louvain'].isin(['5', '1'])]
    assert len(np.intersect1d(ids_summary, matched_data.obs.index)) == matched_data.shape[0], '{}, {}'.format(
        len(ids_summary), matched_data.shape[0])


def test_measure_filter(dataset_api, input_dataset, test_data):
    process_results = handle_selection_ids(dataset_api=dataset_api, dataset=input_dataset,
        data_filter={'filters': [['DSCR3', '>', 2]]})
    ids_summary = process_results['ids']
    matched_data = test_data[test_data[:, 'DSCR3'].X > 2]
    assert len(np.intersect1d(ids_summary, matched_data.obs.index)) == matched_data.shape[0]
