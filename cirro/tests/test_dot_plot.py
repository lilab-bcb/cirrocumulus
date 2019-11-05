import os
import unittest

import anndata
import numpy as np
import pandas as pd

from cirro.data_processing import process_data
from cirro.dataset_api import DatasetAPI
from cirro.entity import Entity
from cirro.parquet_dataset import ParquetDataset

measures = ['DSCR3', 'SUMO3', 'TNFRSF4', 'JAKMIP1']
by = 'louvain'
path = 'test-data/3K_PBMC.pq'
data = anndata.read('test-data/3K_PBMC.h5ad')
dataset = Entity(path, {'name': os.path.splitext(os.path.basename(path))[0], 'url': path})
pq_dataset = ParquetDataset()


class TestDotPlot(unittest.TestCase):

    def setUp(self):
        self.dataset_api = DatasetAPI()
        self.dataset_api.add(['pq'], pq_dataset)

    def test_dot_plot(self):
        X = data[:, measures].X
        X = X.toarray()
        df = pd.DataFrame(data=X, columns=measures)
        df[by] = data.obs[by].values

        def fraction_expressed(g):
            return (g > 0).sum() / len(g)

        summarized_df = df.groupby(by).aggregate(['mean', fraction_expressed])
        process_results = process_data(dataset_api=self.dataset_api, dataset=dataset, dotplot_measures=measures,
            dotplot_dimensions=[by], return_types=['dotplot'])
        dotplot_results = process_results['dotplot'].collect()[0]
        values = dotplot_results['values']
        for key in measures:
            index = -1
            for i in range(len(values)):
                if values[i]['name'] == key:
                    index = i
                    break
            np.testing.assert_allclose(summarized_df[key]['mean'].values, values[index]['mean'],
                atol=0.0000001, err_msg='mean')
            np.testing.assert_allclose(summarized_df[key]['fraction_expressed'].values,
                values[index]['fraction_expressed'], err_msg='fraction_expressed')


if __name__ == "__main__":
    unittest.main()
