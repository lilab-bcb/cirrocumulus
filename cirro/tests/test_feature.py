import os
import unittest

import anndata
import numpy as np

from cirro.data_processing import process_data
from cirro.dataset_api import DatasetAPI
from cirro.entity import Entity
from cirro.parquet_dataset import ParquetDataset

measures = ['DSCR3', 'SUMO3', 'TNFRSF4']
dimensions = ['louvain']
path = 'test-data/3K_PBMC.pq'
data = anndata.read('test-data/3K_PBMC.h5ad')
dataset = Entity(path, {'name': os.path.splitext(os.path.basename(path))[0], 'url': path})
pq_dataset = ParquetDataset()


class TestFeature(unittest.TestCase):

    def setUp(self):
        self.dataset_api = DatasetAPI()
        self.dataset_api.add(['pq'], pq_dataset)

    def diff_measures(self, measures_summary):
        mean = data[:, measures].X.mean(axis=0)
        min = data[:, measures].X.min(axis=0)
        max = data[:, measures].X.max(axis=0)
        fraction_expressed = (data[:, measures].X > 0).sum(axis=0) / data.shape[0]
        np.testing.assert_allclose(
            measures_summary.iloc[:, measures_summary.columns.get_level_values(1) == 'min'].values.flatten(), min,
            err_msg='min')
        np.testing.assert_allclose(
            measures_summary.iloc[:, measures_summary.columns.get_level_values(1) == 'max'].values.flatten(), max,
            err_msg='max')
        np.testing.assert_allclose(
            measures_summary.iloc[:,
            measures_summary.columns.get_level_values(1) == 'fraction_expressed'].values.flatten(),
            fraction_expressed, err_msg='fraction_expressed')
        np.testing.assert_allclose(
            measures_summary.iloc[:, measures_summary.columns.get_level_values(1) == 'mean'].values.flatten(), mean,
            err_msg='mean', atol=0.0000001)

    def diff_dimensions(self, dimension_value_counts):
        for key in dimensions:
            value_counts = data.obs[key].value_counts()
            computed_values_counts = dimension_value_counts[key]
            self.assertEqual((computed_values_counts[computed_values_counts.columns[0]] - value_counts).sum(), 0, key)

    def test_measure_and_dimension(self):
        process_results = process_data(dataset_api=self.dataset_api, dataset=dataset, summary_measures=measures,
            summary_dimensions=dimensions, return_types=['summary'])

        measures_summary, dimension_value_counts = process_results['summary'].collect()
        self.diff_measures(measures_summary)
        self.diff_dimensions(dimension_value_counts)

    def test_measure(self):
        process_results = process_data(dataset_api=self.dataset_api, dataset=dataset, summary_measures=measures,
            return_types=['summary'])
        measures_summary, dimension_value_counts = process_results['summary'].collect()
        self.diff_measures(measures_summary)

    def test_dimension(self):
        process_results = process_data(dataset_api=self.dataset_api, dataset=dataset,
            summary_dimensions=dimensions, return_types=['summary'])
        measures_summary, dimension_value_counts = process_results['summary'].collect()
        self.diff_dimensions(dimension_value_counts)


if __name__ == "__main__":
    unittest.main()
