import anndata
import pandas as pd
import scipy


class H5ADBackend:

    def __init__(self, backed='r'):
        self.path_to_data = {}
        self.backed = backed
        # only works with local files

    def get_df(self, filesystem, path, keys, embedding_key=None, index=False):
        adata = self.path_to_data.get(path, None)
        if adata is None:
            adata = anndata.read(path, backed=self.backed)
            self.path_to_data[path] = adata
        return self.__get_df(adata, keys, embedding_key, index=index)

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
        result['var'] = adata.var_names.values
        result['obs'] = obs
        result['obs_cat'] = obs_cat
        result['n_obs'] = adata.shape[0]
        embeddings = []
        for key in adata.obsm_keys():
            if key.startswith('X_'):
                embeddings.append(dict(name=key, dimensions=adata.obsm[key].shape[1]))
        result['embeddings'] = embeddings
        return result

    def __get_df(self, adata, keys, embedding_key=None, index=False):
        is_obs = True
        if index:
            df = pd.DataFrame(index=(adata.obs.index.values if is_obs else adata.var.index.values))
        else:
            df = pd.DataFrame(index=pd.RangeIndex(adata.shape[0 if is_obs else 1]))
        for i in range(len(keys)):
            key = keys[i]
            if key in adata.var_names and is_obs:
                X = adata.obs_vector(key)
                if scipy.sparse.issparse(X):
                    X = X.toarray()
                values = X
            elif key in adata.obs and is_obs:
                values = adata.obs[key].values
            else:
                raise ValueError('{} not found'.format(key))
            df[key] = values
        if embedding_key is not None:
            embedding_name = embedding_key['name']
            embedding_data = adata.obsm[embedding_name]
            dimensions = embedding_key['dimensions']
            for i in range(dimensions):
                df[embedding_key['coordinate_columns'][i]] = embedding_data[:, i]
        return df
