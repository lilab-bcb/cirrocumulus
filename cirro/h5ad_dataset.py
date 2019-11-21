import anndata
import pandas as pd
import pyarrow as pa
import scipy


class H5ADDataset:

    def __init__(self, backed='r'):
        self.path_to_data = {}
        self.backed = backed
        # only works with local files

    def schema(self, filesystem, path):
        adata = self.path_to_data.get(path, None)
        if adata is None:
            adata = anndata.read(path, backed=self.backed)
            self.path_to_data[path] = adata
        obs_cat = []
        obs = []
        result = {'version': '1'}
        for key in adata.obs_keys():
            if pd.api.types.is_categorical_dtype(adata.obs[key]):
                obs_cat.append(key)
            else:
                obs.append(key)
        result['var'] = adata.var_names.values.tolist()
        result['obs'] = obs
        result['obsCat'] = obs_cat
        result['nObs'] = adata.shape[0]
        embeddings = []
        for key in adata.obsm_keys():
            if key.startswith('X_'):
                embeddings.append(dict(name=key, dimensions=adata.obsm[key].shape[1]))
        result['embeddings'] = embeddings
        return result

    def statistics(self, file_system, path, keys, basis):
        py_dict = self.get_py_dict(file_system, path, keys, basis)
        key_to_stats = {}
        for key in py_dict:
            key_to_stats[key] = (py_dict[key].min(), py_dict[key].max())
        return key_to_stats

    def table(self, file_system, path, keys, basis=None):
        return pa.Table.from_pydict(self.get_py_dict(file_system, path, keys, basis=basis))

    def tables(self, file_system, path, keys, basis=None):
        yield self.table(file_system, path, keys, basis)

    def get_py_dict(self, file_system, path, keys, basis=None):
        is_obs = True
        py_dict = {}
        adata = self.path_to_data.get(path, None)
        if adata is None:
            adata = anndata.read(path, backed=self.backed)
            self.path_to_data[path] = adata

        for i in range(len(keys)):
            key = keys[i]
            values = None
            if key == 'index':
                values = adata.obs.index.values if is_obs else adata.var.index.values
            else:
                if key in adata.var_names and is_obs:
                    X = adata.obs_vector(key)
                    if scipy.sparse.issparse(X):
                        X = X.toarray()
                    values = X
                elif key in adata.obs and is_obs:
                    values = adata.obs[key].values
                else:
                    print('{} not found'.format(key))
            if values is not None:
                py_dict[key] = values
        if basis is not None:
            embedding_name = basis['name']
            embedding_data = adata.obsm[embedding_name]
            dimensions = basis['dimensions']
            for i in range(dimensions):
                py_dict[basis['coordinate_columns'][i]] = embedding_data[:, i]
        return py_dict
