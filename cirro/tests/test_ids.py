import os
import unittest

import anndata
import numpy as np
import pandas as pd

from cirro.data_processing import process_data
from cirro.dataset_api import DatasetAPI
from cirro.embedding_aggregator import get_basis, EmbeddingAggregator
from cirro.entity import Entity
from cirro.h5ad_dataset import H5ADDataset

nbins = 100
measures = ['DSCR3', 'SUMO3', 'TNFRSF4']
dimensions = ['louvain']
path = 'test-data/3K_PBMC.h5ad'
basis = get_basis('X_umap')
dataset = Entity(path, {'id': path, 'name': os.path.splitext(os.path.basename(path))[0], 'url': path})


class TestIds(unittest.TestCase):

    def setUp(self):
        self.dataset_api = DatasetAPI()
        self.dataset_api.add(H5ADDataset(backed=None))
        data = anndata.read(path)
        data.obs = data.obs.reset_index()
        self.data = data

    def test_dimension_filter(self):
        process_results = process_data(dataset_api=self.dataset_api, dataset=dataset, return_types=['ids'],
            data_filter={'filters': [['obs.louvain', 'in', ['1', '5']]]})
        ids_summary = process_results['ids']
        matched_data = self.data[self.data.obs['louvain'].isin(['5', '1'])]
        self.assertEqual(len(np.intersect1d(ids_summary, matched_data.obs['index'])), matched_data.shape[0],
            '{}, {}'.format(len(ids_summary), matched_data.shape[0]))

    def test_binning_dimension_filter(self):
        df = pd.DataFrame(self.data[:, measures].X.toarray(), columns=measures)
        df = df.join(self.data.obs[dimensions])
        df['index'] = self.data.obs['index']
        df = df.join(pd.DataFrame(self.data.obsm['X_umap'][:, 0:2], columns=basis['coordinate_columns']))
        EmbeddingAggregator.convert_coords_to_bin(df, nbins=nbins,
            coordinate_columns=basis['coordinate_columns'],
            coordinate_column_to_range=None)
        df = df[df['louvain'].isin(['5', '1'])]
        process_results = process_data(dataset_api=self.dataset_api, basis=basis, nbins=nbins, dataset=dataset,
            return_types=['ids'], data_filter={'filters': [['obs.louvain', 'in', ['1', '5']]]})
        ids_summary = process_results['ids']
        both = np.intersect1d(ids_summary, df['index'])
        self.assertEqual(len(both), df.shape[0], '{}, {}'.format(len(ids_summary), len(df)))

    def test_measure_filter(self):
        process_results = process_data(dataset_api=self.dataset_api, dataset=dataset, return_types=['ids'],
            data_filter={'filters': [['DSCR3', '>', 2]]})
        ids_summary = process_results['ids']
        matched_data = self.data[self.data[:, 'DSCR3'].X > 2]
        self.assertEqual(len(np.intersect1d(ids_summary, matched_data.obs['index'])), matched_data.shape[0])

    def test_binning_measure_filter(self):
        df = pd.DataFrame(self.data[:, measures].X.toarray(), columns=measures)
        df = df.join(self.data.obs[dimensions])
        df['index'] = self.data.obs['index']
        df = df.join(pd.DataFrame(self.data.obsm['X_umap'][:, 0:2], columns=basis['coordinate_columns']))
        EmbeddingAggregator.convert_coords_to_bin(df, nbins=nbins,
            coordinate_columns=basis['coordinate_columns'],
            coordinate_column_to_range=None)
        df = df[df['DSCR3'] > 2]
        process_results = process_data(dataset_api=self.dataset_api, basis=basis, nbins=nbins, dataset=dataset,
            return_types=['ids'], data_filter={'filters': [['DSCR3', '>', 2]]})
        ids_summary = process_results['ids']
        both = np.intersect1d(ids_summary, df['index'])
        self.assertEqual(len(both), df.shape[0])


if __name__ == "__main__":
    unittest.main()
