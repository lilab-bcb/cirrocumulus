import numpy as np
import pandas as pd
import scipy.sparse


class SimpleData:

    def __init__(self, X, obs, var):
        self.obs = obs
        self.var = var
        self.uns = {}
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
        if len(df) == 0:
            zeros = np.full(len(var_ids), 0)
            empty = np.full(len(var_ids), np.nan)
            return pd.DataFrame(
                data={'min': empty, 'max': empty, 'sum': zeros,
                      'numExpressed': zeros,
                      'mean': empty}, index=var_ids)

        return pd.DataFrame(
            data={'min': np.min(df.values, axis=0), 'max': np.max(df.values, axis=0), 'sum': df.sum().values,
                  'numExpressed': (df.values != 0).sum(axis=0),
                  'mean': df.mean().values}, index=var_ids)

    @staticmethod
    def get_var_indices(adata, names):
        return adata.var.index.get_indexer_for(names)

    @staticmethod
    def find_markers(adata, key, n_genes):
        import scipy.stats as ss
        marker_results = []  # array of category, name, id, features
        category = key + ' markers'

        for cat in adata.obs[key].cat.categories:
            print('Computing markers for {}, {}'.format(key, cat))
            mask = adata.obs[key] == cat
            ds1 = adata[mask]
            ds_rest = adata[~mask]
            gene_names_keep = ds1.X.mean(axis=0) > ds_rest.X.mean(axis=0)
            if isinstance(gene_names_keep, np.matrix):
                gene_names_keep = gene_names_keep.A1
            ds1 = ds1[:, gene_names_keep]
            ds_rest = ds_rest[:, gene_names_keep]
            pvals = np.full(ds1.shape[1], 1.0)
            ds1_X = ds1.X
            ds_rest_X = ds_rest.X

            for i in range(ds1.shape[1]):
                v1 = ds1_X[:, i]
                v2 = ds_rest_X[:, i]
                if v1.data.size > 0 and v2.data.size > 0:
                    try:
                        _, pvals[i] = ss.mannwhitneyu(v1.toarray()[:, 0], v2.toarray()[:, 0],
                            alternative="two-sided")
                    except ValueError:
                        # All numbers are identical
                        pass
            fc = ds1_X.mean(axis=0) - ds_rest_X.mean(axis=0)
            if isinstance(fc, np.matrix):
                fc = fc.A1
            df = pd.DataFrame(data=dict(pvals=pvals, fc=fc), index=ds1.var_names)
            df = df.sort_values(by=['pvals', 'fc'], ascending=[True, False])
            features = df[:n_genes].index.values
            marker_results.append(dict(category=category, name=str(cat), features=features))
        return marker_results

    @staticmethod
    def has_markers(adata):
        return hasattr(adata, 'uns') and 'rank_genes_groups' in adata.uns or 'de_res' in adata.varm

    @staticmethod
    def schema(adata):
        obs_cat = []
        obs = []
        result = {'version': '1.0.0'}
        marker_results = []
        marker_results += adata.uns.get('markers', [])
        result['markers'] = marker_results
        n_genes = 10
        if SimpleData.has_markers(adata):
            if hasattr(adata, 'uns') and 'rank_genes_groups' in adata.uns:  # scanpy
                for key in adata.uns.keys():
                    rank_genes_groups = adata.uns[key]
                    if isinstance(rank_genes_groups, dict) and 'logfoldchanges' in rank_genes_groups:
                        groupby = str(rank_genes_groups['params']['groupby'])
                        group_names = rank_genes_groups['names'].dtype.names
                        category = groupby + ' markers'
                        for group_name in group_names:
                            gene_names = rank_genes_groups['names'][group_name]
                            # scores = rank_genes_groups['scores'][group_name]
                            features = gene_names[:n_genes]
                            marker_results.append(dict(category=category, name=str(group_name), features=features))
            else:  # pegasus
                de_res = adata.varm['de_res']
                names = de_res.dtype.names
                field_names = set()  # e.g. 1:auroc
                cluster_names = set()
                for name in names:
                    index = name.index(':')
                    field_name = name[index + 1:]
                    cluster_name = name[:index]
                    field_names.add(field_name)
                    cluster_names.add(cluster_name)
                sort_field_names = ['mwu_qval', 'auroc', 't_qval']
                de_res = pd.DataFrame(data=de_res, index=adata.var.index)
                field_use = None
                for field in sort_field_names:
                    if field in field_names:
                        field_use = field
                        break
                if field_use is not None:
                    field_ascending = field_use != 'auroc'
                    for cluster_name in cluster_names:
                        fc_column = '{}:log2FC'.format(cluster_name)
                        name = '{}:{}'.format(cluster_name, field_name)
                        idx_up = de_res[fc_column].values > 0
                        df_up = de_res.loc[idx_up].sort_values(by=[name, fc_column], ascending=[field_ascending, False])
                        features = df_up[:n_genes].index.values
                        marker_results.append(dict(category='markers', name=str(cluster_name), features=features))

        for key in adata.obs_keys():
            if pd.api.types.is_categorical_dtype(adata.obs[key]) or pd.api.types.is_bool_dtype(
                    adata.obs[key]) or pd.api.types.is_object_dtype(adata.obs[key]):
                obs_cat.append(key)
            else:
                obs.append(key)
        # spatial_node = adata.uns['spatial'] if 'spatial' in adata.uns else None
        #
        # if spatial_node is not None:
        #     spatial_node_keys = list(spatial_node.keys())  # list of datasets
        #     if len(spatial_node_keys) == 1:
        #         spatial_node = spatial_node[spatial_node_keys[0]]

        images_node = adata.uns.get('images', [])
        image_names = list(map(lambda x: x['name'], images_node))
        result['var'] = adata.var_names.values
        result['obs'] = obs
        result['obsCat'] = obs_cat
        result['shape'] = adata.shape

        embeddings = []
        for key in adata.obsm_keys():
            dim = min(3, adata.obsm[key].shape[1])
            if dim < 2:
                continue
            embedding = dict(name=key, dimensions=dim)
            try:
                image_index = image_names.index(key)
                embedding['spatial'] = images_node[image_index]
            except ValueError:
                pass

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
