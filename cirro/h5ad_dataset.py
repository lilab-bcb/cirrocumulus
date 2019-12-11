import anndata
import pandas as pd
from cirro.simple_data import SimpleData


class H5ADDataset:

    def __init__(self, backed='r'):
        self.path_to_data = {}
        self.backed = backed
        # only works with local files

    def get_suffixes(self):
        return ['h5ad']

    def get_data(self, path):
        adata = self.path_to_data.get(path)
        if adata is None:
            adata = anndata.read(path, backed=self.backed)
            self.path_to_data[path] = adata
        return adata

    def schema(self, filesystem, path):
        return SimpleData.schema(self.get_data(path))

    def statistics(self, file_system, path, keys, basis):
        py_dict = self.get_py_dict(file_system, path, keys, basis)
        key_to_stats = {}
        for key in py_dict:
            key_to_stats[key] = (py_dict[key].min(), py_dict[key].max())
        return key_to_stats

    def read(self, file_system, path, obs_keys=[], var_keys=[], basis=None, dataset=None):
        adata = self.get_data(path)
        obs = None
        X = None

        if len(var_keys) > 0:
            X = adata[:, var_keys].X
        for key in obs_keys:
            if key == 'id':
                values = adata.obs.index.values
            else:
                values = adata.obs[key].values
            if obs is None:
                obs = pd.DataFrame()
            obs[key] = values

        if basis is not None:
            if obs is None:
                obs = pd.DataFrame()
            embedding_name = basis['name']
            embedding_data = adata.obsm[embedding_name]
            dimensions = basis['dimensions']
            for i in range(dimensions):
                obs[basis['coordinate_columns'][i]] = embedding_data[:, i]
        return SimpleData(X, obs, pd.DataFrame(index=pd.Index(var_keys)))
