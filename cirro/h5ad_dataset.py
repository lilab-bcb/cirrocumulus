import anndata
import numpy as np
import pandas as pd
import scipy.sparse

from cirro.simple_data import SimpleData


class H5ADDataset:

    def __init__(self, backed='r', force_sparse=True):
        self.path_to_data = {}
        self.backed = backed
        self.force_sparse = force_sparse
        # only works with local files

    def get_suffixes(self):
        return ['h5ad', 'loom']

    def read_adata(self, path):
        path_lc = path.lower()
        if path_lc.endswith('.loom'):
            return anndata.read_loom(path)
        return anndata.read(path, backed=self.backed)
        # elif path.endswith('.mtx'):
        #
        #     return anndata.read_mtx(path, backed=self.backed)
        # elif path.endswith('.txt'):
        #
        #     return anndata.read_text(path, backed=self.backed)
        # elif path.endswith('.csv'):
        #
        #     return anndata.read_csv(path, backed=self.backed)

    def get_data(self, path):
        adata = self.path_to_data.get(path)
        if adata is None:
            adata = self.read_adata(path)
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
            indexer = adata.var.index.get_indexer_for(var_keys)
            indexer_sort = np.argsort(indexer)
            X = adata[:, indexer[indexer_sort]].X
            if len(X.shape) == 1:
                X = np.array([X]).T

            var_keys = np.array(var_keys)[indexer_sort]
            if self.force_sparse and not scipy.sparse.issparse(X):
                X = scipy.sparse.csr_matrix(X)
        for key in obs_keys:
            if key == 'index':
                values = adata.obs.index.values
            else:
                values = adata.obs[key].values
            if obs is None:
                obs = pd.DataFrame()
            obs[key] = values

        if basis is not None:
            if obs is None:
                obs = pd.DataFrame()
            for b in basis:
                embedding_name = b['name']
                embedding_data = adata.obsm[embedding_name]
                dimensions = b['dimensions']
                for i in range(dimensions):
                    obs[b['coordinate_columns'][i]] = embedding_data[:, i]
        return SimpleData(X, obs, pd.DataFrame(index=pd.Index(var_keys)))
