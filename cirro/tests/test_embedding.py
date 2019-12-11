import os
import unittest

import anndata
import numpy as np
import pandas as pd

from cirro.data_processing import process_data
from cirro.dataset_api import DatasetAPI
from cirro.embedding_aggregator import get_basis, EmbeddingAggregator
from cirro.entity import Entity
from cirro.parquet_dataset import ParquetDataset


measures = ['DSCR3', 'SUMO3', 'TNFRSF4']
dimensions = ['louvain']
path = 'test-data/3K_PBMC.pq'
data = anndata.read('test-data/3K_PBMC.h5ad')
basis = get_basis('X_umap')
dataset = Entity(path, {'name': os.path.splitext(os.path.basename(path))[0], 'url': path})
pq_dataset = ParquetDataset()


class TestEmbedding(unittest.TestCase):

    def setUp(self):
        self.dataset_api = DatasetAPI()
        self.dataset_api.add(pq_dataset)
        df = pd.DataFrame(data[:, measures].X.toarray(), columns=measures)
        df = df.join(data.obs[dimensions].reset_index())
        df = df.join(pd.DataFrame(data.obsm['X_umap'][:, 0:2], columns=basis['coordinate_columns']))
        self.df = df

    def test_no_binning(self):
        process_results = process_data(dataset_api=self.dataset_api, dataset=dataset, basis=basis,
            embedding_measures=measures, embedding_dimensions=dimensions,
            return_types=['embedding'])

        results = process_results['embedding']

        for key in basis['coordinate_columns']:
            np.testing.assert_array_equal(results['coordinates'][key], self.df[key].values, err_msg=key)
        for key in measures:
            np.testing.assert_array_equal(results['values'][key], self.df[key].values, err_msg=key)
        for key in dimensions:
            np.testing.assert_array_equal(results['values'][key], data.obs[key].values, err_msg=key)


    def test_binning(self):
        nbins = 100
        agg_function = 'min'

        agg_dict = {}

        for key in basis['coordinate_columns']:
            agg_dict[key] = 'max'
        for key in measures:
            agg_dict[key] = agg_function
        for key in dimensions:
            agg_dict[key] = lambda x: x.mode()[0]
        df = self.df
        EmbeddingAggregator.convert_coords_to_bin(df, nbins=nbins,
            coordinate_columns=basis['coordinate_columns'],
            coordinate_column_to_range=None)
        # group by bin then agg
        df = df.groupby(df.index).agg(agg_dict)

        process_results = process_data(dataset_api=self.dataset_api, dataset=dataset, basis=basis, nbins=nbins,
            embedding_measures=measures, embedding_dimensions=dimensions,
            return_types=['embedding'], agg_function=agg_function)

        results = process_results['embedding']

        # for key in basis['coordinate_columns']:
        #     np.testing.assert_allclose(measure_df[key].values, df[key].values, atol=0.0000001, err_msg=key)
        # for key in measures:
        #     np.testing.assert_allclose(measure_df[key].values, df[key].values, atol=0.0000001, err_msg=key)
        # for key in dimensions:
        #     np.testing.assert_array_equal(df[key].values, dimension_df[key].values, err_msg=key)
        for key in basis['coordinate_columns']:
            np.testing.assert_array_equal(results['coordinates'][key], df[key].values, err_msg=key)
        for key in measures:
            np.testing.assert_array_equal(results['values'][key], df[key].values, err_msg=key)
        for key in dimensions:
            np.testing.assert_array_equal(results['values'][key]['mode'], df[key].values, err_msg=key)


if __name__ == "__main__":
    unittest.main()
