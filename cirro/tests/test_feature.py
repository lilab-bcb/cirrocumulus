import os
import unittest

import anndata
import pandas as pd

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
        self.dataset_api.add(pq_dataset)

    def diff_measures(self, summary, fields):
        for field in fields:
            my_mean = data[:, field].X.mean()
            my_min = data[:, field].X.min()
            my_max = data[:, field].X.max()
            my_expressed = (data[:, field].X > 0).sum()
            field_summary = summary[field]
            self.assertAlmostEqual(my_mean, field_summary['mean'])
            self.assertAlmostEqual(my_min, field_summary['min'])
            self.assertAlmostEqual(my_max, field_summary['max'])
            self.assertEqual(my_expressed, field_summary['numExpressed'])

    def diff_dimensions(self, summary, fields):
        for key in fields:
            value_counts = data.obs[key].value_counts()
            dimension_summary = summary[key]
            df = pd.DataFrame.from_dict(dimension_summary)
            df['categories'] = df['categories'].astype('category')
            df = df.set_index('categories')
            df = df.iloc[df.index.get_indexer_for(value_counts.index)]
            self.assertEqual((df['counts'] - value_counts).sum(), 0, key)

    def test_measure_and_dimension(self):
        process_results = process_data(dataset_api=self.dataset_api, dataset=dataset, summary_measures=measures,
            summary_dimensions=dimensions, return_types=['summary'])

        summary = process_results['summary']
        self.diff_measures(summary, measures)
        self.diff_dimensions(summary, dimensions)

    def test_measure(self):
        process_results = process_data(dataset_api=self.dataset_api, dataset=dataset, summary_measures=measures,
            return_types=['summary'])
        summary = process_results['summary']
        self.diff_measures(summary, measures)

    def test_dimension(self):
        process_results = process_data(dataset_api=self.dataset_api, dataset=dataset,
            summary_dimensions=dimensions, return_types=['summary'])
        summary = process_results['summary']
        self.diff_dimensions(summary, dimensions)


if __name__ == "__main__":
    unittest.main()
