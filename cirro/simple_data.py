import pandas as pd

import scipy.sparse


class SimpleData:

    def __init__(self, X, obs, var):
        self.X = X
        self.obs = obs
        self.var = var
        n_obs = X.shape[0] if X is not None else (len(obs) if obs is not None else 0)
        n_var = X.shape[1] if X is not None else (len(var) if var is not None else 0)
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
        return df.agg(['min', 'max', 'sum']).T

    @staticmethod
    def X_stats(adata, var_ids):
        indices = SimpleData.get_var_indices(adata, var_ids)
        X = adata.X[:, indices]
        min_values = X.min(axis=0)
        max_values = X.max(axis=0)
        sums = X.sum(axis=0)
        num_expressed = (X > 0).sum(axis=0)
        if scipy.sparse.issparse(X):
            min_values = min_values.toarray().flatten()
            max_values = max_values.toarray().flatten()
            sums = sums.A1
            num_expressed = num_expressed.A1

        return pd.DataFrame(data={'min': min_values, 'max': max_values, 'sum': sums, 'numExpressed': num_expressed},
            index=var_ids)


    # @staticmethod
    # def filter(adata, filter_expr):
    #     return adata.var_index.get_indexer_for(names)

    @staticmethod
    def get_var_indices(adata, names):
        return adata.var.index.get_indexer_for(names)

    @staticmethod
    def get_var_index(adata, name):
        return adata.var.index.get_loc(name)

    @staticmethod
    def schema(adata):
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
            if key.startswith('X_') and adata.obsm[key].shape[1] <= 10:
                embeddings.append(dict(name=key, dimensions=adata.obsm[key].shape[1]))
        result['embeddings'] = embeddings
        return result

    @staticmethod
    def to_df(adata, obs_measures, var_measures, dimensions, basis=None):
        df = pd.DataFrame()
        for key in obs_measures + dimensions:
            df[key] = adata.obs[key]
        indices = adata.get_var_indices(var_measures)
        for i in range(len(var_measures)):
            X = adata.X[:, indices[i]]
            if scipy.sparse.issparse(X):
                X = X.toarray()
            df[key] = X
        if basis is not None:
            embedding_name = basis['name']
            embedding_data = adata.obsm[embedding_name]
            dimensions = basis['dimensions']
            for i in range(dimensions):
                df[basis['coordinate_columns'][i]] = embedding_data[:, i]
        return df
