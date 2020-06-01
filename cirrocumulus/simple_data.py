import numpy as np
import pandas as pd
import scipy.sparse


class SimpleData:

    def __init__(self, X, obs, var):

        self.obs = obs
        self.var = var
        if X is not None:
            if len(X.shape) == 1:
                X = np.array([X]).T
            n_var = X.shape[1]
            n_obs = X.shape[0]
        else:
            n_var = len(var) if var is not None else 0
            n_obs = len(obs) if obs is not None else 0
        self.X = X
        self.shape = (n_obs, n_var)

    @staticmethod
    def view(adata, row_slice):
        X = adata.X[row_slice] if adata.X is not None else None
        obs = adata.obs[row_slice] if adata.obs is not None else None
        return SimpleData(X, obs, adata.var)


    @staticmethod
    def obs_stats(adata, columns):
        df = adata.obs[columns]
        # variables on columns, stats on rows, transpose so that stats are on columns
        return df.agg(['min', 'max', 'sum', 'mean']).T

    @staticmethod
    def X_stats(adata, var_ids):
        indices = SimpleData.get_var_indices(adata, var_ids)

        X = adata.X[:, indices]
        min_values = X.min(axis=0)
        mean_values = X.mean(axis=0)
        max_values = X.max(axis=0)
        sums = X.sum(axis=0)

        if scipy.sparse.issparse(X):
            min_values = min_values.toarray().flatten()
            max_values = max_values.toarray().flatten()
            mean_values = mean_values.A1
            sums = sums.A1
            num_expressed = X.getnnz(axis=0)
        else:
            num_expressed = (X != 0).sum(axis=0)

        return pd.DataFrame(data={'min': min_values, 'max': max_values, 'sum': sums, 'numExpressed': num_expressed,
                                  'mean': mean_values}, index=var_ids)

    @staticmethod
    def get_var_indices(adata, names):
        return adata.var.index.get_indexer_for(names)

    @staticmethod
    def schema(adata):
        obs_cat = []
        obs = []
        result = {'version': '1.0.0'}
        for key in adata.obs_keys():
            if pd.api.types.is_categorical_dtype(adata.obs[key]) or pd.api.types.is_bool_dtype(adata.obs[key]):
                obs_cat.append(key)
            else:
                obs.append(key)
        spatial_node = adata.uns['spatial'] if 'spatial' in adata.uns else None

        if spatial_node is not None:
            spatial_node_keys = list(spatial_node.keys())
            if len(spatial_node_keys) == 1:
                spatial_node = spatial_node[spatial_node_keys[0]].keys()  # images', 'metadata', 'scalefactors']

        result['var'] = list(
            sorted(adata.var_names.values, key=lambda x: ('zzzzz' + x.lower()) if x[0].isdigit() else x.lower()))
        result['obs'] = obs
        result['obsCat'] = obs_cat
        result['shape'] = adata.shape
        embeddings = []
        for key in adata.obsm_keys():

            dim = min(3, adata.obsm[key].shape[1])
            embedding = dict(name=key, dimensions=dim)

            if key == 'spatial':
                if spatial_node is not None and 'scalefactors' in spatial_node and 'images' in spatial_node:
                    embedding['spatial'] = dict(scalefactors=spatial_node['scalefactors'],
                        images=list(spatial_node['images'].keys()))
            if dim == 3:
                embeddings.append(embedding)
                embedding = embedding.copy()
                embedding['dimensions'] = 2
                embeddings.append(embedding)
            else:
                embeddings.append(embedding)
        result['embeddings'] = embeddings
        return result

    @staticmethod
    def to_df(adata, obs_measures, var_measures, dimensions, basis=None):
        df = pd.DataFrame()
        obs_keys = obs_measures + dimensions
        if basis is not None:
            obs_keys += basis['coordinate_columns']

        for key in obs_keys:
            df[key] = adata.obs[key]
        indices = SimpleData.get_var_indices(adata, var_measures)
        for i in range(len(var_measures)):
            X = adata.X[:, indices[i]]
            if scipy.sparse.issparse(X):
                X = X.toarray()
            df[var_measures[i]] = X
        return df
