import anndata
import numpy as np
import pandas as pd
import scipy.sparse

from cirrocumulus.abstract_dataset import AbstractDataset
from cirrocumulus.simple_data import SimpleData


class AnndataDataset(AbstractDataset):

    def __init__(self, backed=None, force_sparse=True, extensions=['h5ad', 'loom', 'zarr']):
        super().__init__()
        self.path_to_data = {}
        self.backed = backed
        self.force_sparse = force_sparse
        self.extensions = extensions
        self.cached_dataset_id = None
        self.cached_data = {}
        # only works with local files

    def get_suffixes(self):
        return self.extensions

    def read_adata(self, path):
        path_lc = path.lower()
        if path_lc.endswith('.loom'):
            return anndata.read_loom(path)
        elif path_lc.endswith('.zarr'):
            return anndata.read_zarr(path)
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

    def add_data(self, path, data):
        self.path_to_data[path] = data

    def get_data(self, path):
        adata = self.path_to_data.get(path)
        if adata is None:
            adata = self.read_adata(path)
            if scipy.sparse.isspmatrix_csr(adata.X):
                adata.X = adata.X.tocsc()
            self.add_data(path, adata)
        return adata

    def schema(self, fs_adapter, path):
        return SimpleData.schema(self.get_data(path))

    def read(self, fs_adapter, path, obs_keys=[], var_keys=[], basis=None, dataset=None, schema=None):
        adata = self.get_data(path)
        obs = None
        X = None

        if self.cached_dataset_id != dataset.id:
            self.cached_dataset_id = dataset.id
            self.cached_data = {}

        if len(var_keys) > 0:
            cache_var_key = str(dataset.id) + '-var_keys'
            cache_X_key = str(dataset.id) + '-X'
            cached_var_keys = self.cached_data.get(cache_var_key)

            if var_keys == cached_var_keys:
                X = self.cached_data[cache_X_key]
            else:
                X = adata[:, var_keys].X
                if len(X.shape) == 1:
                    X = np.array([X]).T
                if self.force_sparse and not scipy.sparse.issparse(X):
                    X = scipy.sparse.csr_matrix(X)

                self.cached_data[cache_X_key] = X
                self.cached_data[cache_var_key] = var_keys

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
