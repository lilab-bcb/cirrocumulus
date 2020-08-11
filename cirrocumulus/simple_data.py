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
    def obs_stats(df, columns):
        df = df[columns]
        # variables on columns, stats on rows, transpose so that stats are on columns
        return df.agg(['min', 'max', 'sum', 'mean']).T

    @staticmethod
    def X_stats(df, var_ids):
        df = df[var_ids]
        return pd.DataFrame(
            data={'min': np.min(df.values, axis=0), 'max': np.max(df.values, axis=0), 'sum': df.sum().values,
                  'numExpressed': (df.values != 0).sum(axis=0),
                  'mean': df.mean().values}, index=var_ids)


    @staticmethod
    def get_var_indices(adata, names):
        return adata.var.index.get_indexer_for(names)

    @staticmethod
    def find_markers(adata, key, marker_dict, n_genes):
        import scipy.stats as ss
        nfeatures = adata.shape[1]
        gene_names = adata.var_names
        markers = {}
        marker_dict[key + ' markers'] = markers
        for cat in adata.obs[key].cat.categories:
            stats = np.zeros(nfeatures, dtype=np.float32)
            pvals = np.full(nfeatures, 1.0)
            mask = adata.obs[key] == cat
            ds1 = adata[mask]
            ds_rest = adata[~mask]
            for i in range(nfeatures):
                v1 = ds1.X[:, i]
                v2 = ds_rest.X[:, i]
                if v1.data.size > 0 and v2.data.size > 0:
                    stats[i], pvals[i] = ss.mannwhitneyu(v1.toarray()[:, 0], v2.toarray()[:, 0],
                        alternative="two-sided")
            keep = stats > 0
            pvals = pvals[keep]
            order = np.argsort(pvals)
            markers[str(cat)] = gene_names[order][:n_genes]

    @staticmethod
    def has_markers(adata):
        return hasattr(adata, 'uns') and 'rank_genes_groups' in adata.uns

    @staticmethod
    def schema(adata):
        obs_cat = []
        obs = []
        result = {'version': '1.0.0'}
        marker_dict = adata.uns.get('markers', {})
        result['markers'] = marker_dict
        if 'seurat_clusters' in adata.obs:
            adata.obs['seurat_clusters'] = adata.obs['seurat_clusters'].astype('category')
        if SimpleData.has_markers(adata):
            rank_genes_groups = adata.uns['rank_genes_groups']
            groupby = str(rank_genes_groups['params']['groupby'])
            group_names = rank_genes_groups['names'].dtype.names
            n_genes = 10
            markers = {}
            marker_dict[groupby + ' markers'] = markers
            for group_name in group_names:
                gene_names = rank_genes_groups['names'][group_name]
                # scores = rank_genes_groups['scores'][group_name]
                markers[group_name] = gene_names[:n_genes]
        for key in adata.obs_keys():
            if pd.api.types.is_categorical_dtype(adata.obs[key]) or pd.api.types.is_bool_dtype(
                    adata.obs[key]) or pd.api.types.is_object_dtype(adata.obs[key]):
                obs_cat.append(key)
            else:
                obs.append(key)
        spatial_node = adata.uns['spatial'] if 'spatial' in adata.uns else None

        if spatial_node is not None:
            spatial_node_keys = list(spatial_node.keys())
            if len(spatial_node_keys) == 1:
                spatial_node = spatial_node[spatial_node_keys[0]].keys()  # images', 'metadata', 'scalefactors']

        result['var'] = adata.var_names.values
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
