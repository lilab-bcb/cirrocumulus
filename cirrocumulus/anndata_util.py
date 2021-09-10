from enum import Enum

import numpy as np
import pandas as pd
from cirrocumulus.io_util import cirro_id


class DataType(Enum):
    UNSPECIFIED = 'unspecified'
    MODULE = 'module'


def get_scanpy_marker_keys(dataset):
    scanpy_marker_keys = []
    for key in dataset.uns.keys():
        rank_genes_groups = dataset.uns[key]
        if isinstance(rank_genes_groups, dict) and 'names' in rank_genes_groups and (
                'pvals' in rank_genes_groups or 'pvals_adj' in rank_genes_groups or 'scores' in rank_genes_groups):
            scanpy_marker_keys.append(key)
    return scanpy_marker_keys


def get_pegasus_marker_keys(dataset):
    marker_keys = []
    for key in dataset.varm.keys():
        d = dataset.varm[key]
        if isinstance(d, np.recarray):
            try:
                ''.join(d.dtype.names).index('log2Mean')
                marker_keys.append(key)
            except ValueError:
                pass
    return marker_keys


def obs_stats(adata, columns):
    df = adata.obs[columns]
    # variables on columns, stats on rows, transpose so that stats are on columns
    return df.agg(['min', 'max', 'sum', 'mean']).T


def X_stats(adata):
    X = adata.X
    return pd.DataFrame(
        data={'min': X.min(axis=0).toarray().flatten(),
              'max': X.max(axis=0).toarray().flatten(),
              'sum': X.sum(axis=0).flatten(),
              'numExpressed': X.getnnz(axis=0),
              'mean': X.mean(axis=0)}, index=adata.var.index)


def datasets_schema(datasets):
    """ Gets dataset schema.

    Returns
        schema dict. Example:
        {"version":"1.0.0",
        "categoryOrder":{
            "louvain":["0","1","2","3","4","5","6","7"],
            "leiden":["0","1","2","3","4","5","6","7"]},
        "var":["TNFRSF4","CPSF3L","ATAD3C"],
        "obs":["percent_mito","n_counts"],
        "obsCat":["louvain","leiden"],
        "shape":[2638,1838],
        "embeddings":[{"name":"X_pca","dimensions":3},{"name":"X_pca","dimensions":2},{"name":"X_umap","dimensions":2}]
    }
    """
    obs_cat = []
    obs = []

    marker_results = []
    de_results_format = 'records'
    for dataset in datasets:
        prior_marker_results = dataset.uns.get('markers', [])

        if isinstance(prior_marker_results, str):
            import json
            prior_marker_results = json.loads(prior_marker_results)
    marker_results += prior_marker_results

    n_genes = 10
    de_results = []  # array of dicts containing params logfoldchanges, pvals_adj, scores, names
    category_to_order = {}
    embeddings = []
    field_to_value_to_color = dict()  # field -> value -> color
    for i in range(len(datasets)):
        dataset = datasets[i]
        for scanpy_marker_key in get_scanpy_marker_keys(dataset):
            rank_genes_groups = dataset.uns[scanpy_marker_key]
            has_fc = 'logfoldchanges' in rank_genes_groups
            min_fold_change = 1
            params = rank_genes_groups['params']
            prefix = None if i == 0 else dataset.uns.get('name')

            # pts and pts_rest in later scanpy versions
            rank_genes_groups_keys = list(rank_genes_groups.keys())
            for k in ['params', 'names']:
                if k in rank_genes_groups_keys:
                    rank_genes_groups_keys.remove(k)
            if 'pvals' in rank_genes_groups_keys and 'pvals_adj' in rank_genes_groups_keys:
                rank_genes_groups_keys.remove('pvals')
            category = '{} ({})'.format(params['groupby'], scanpy_marker_key)
            de_result_name = category if prefix is None else prefix + '-' + category
            de_result_df = None
            group_names = rank_genes_groups['names'].dtype.names
            for group_name in group_names:
                group_df = pd.DataFrame(index=rank_genes_groups['names'][group_name])
                group_df = group_df[group_df.index != 'nan']
                for rank_genes_groups_key in rank_genes_groups_keys:
                    values = rank_genes_groups[rank_genes_groups_key][group_name]
                    column_name = '{}:{}'.format(group_name, rank_genes_groups_key)
                    group_df[column_name] = values

                if de_result_df is None:
                    de_result_df = group_df
                else:
                    de_result_df = de_result_df.join(group_df, how='outer')
                if n_genes > 0:
                    markers_df = group_df
                    if has_fc:
                        markers_df = group_df[group_df['{}:logfoldchanges'.format(group_name)] > min_fold_change]
                    if len(markers_df) > 0:
                        marker_results.append(
                            dict(category=de_result_name, name=str(group_name), features=markers_df.index[:n_genes]))

            if de_results_format == 'records':
                de_result_data = de_result_df.reset_index().to_dict(orient='records')
            else:
                de_result_data = dict(index=de_result_df.index)
                for c in de_result_df:
                    de_result_data[c] = de_result_df[c]

            de_result = dict(id=cirro_id(),
                             color='logfoldchanges' if 'logfoldchanges' in rank_genes_groups_keys else
                             rank_genes_groups_keys[0],
                             size='pvals_adj' if 'pvals_adj' in rank_genes_groups_keys else rank_genes_groups_keys[
                                 0],
                             params=params, groups=group_names, fields=rank_genes_groups_keys, type='de',
                             name=de_result_name)
            de_result['data'] = de_result_data
            de_results.append(de_result)
        for varm_field in get_pegasus_marker_keys(dataset):
            de_res = dataset.varm[varm_field]
            names = de_res.dtype.names
            field_names = set()  # e.g. 1:auroc
            group_names = set()
            for name in names:
                index = name.rindex(':')
                field_name = name[index + 1:]
                group_name = name[:index]
                field_names.add(field_name)
                group_names.add(group_name)
            group_names = list(group_names)
            field_names = list(field_names)

            de_result_df = pd.DataFrame(data=de_res, index=dataset.var.index)
            de_result_df.index.name = 'index'
            if de_results_format == 'records':
                de_result_data = de_result_df.reset_index().to_dict(orient='records')
            else:
                de_result_data = dict(index=de_result_df.index)
                for c in de_res:
                    de_result_data[c] = de_result_df[c]

            de_result = dict(id=cirro_id(), type='de', name='pegasus_de',
                             color='log2FC' if 'log2FC' in field_names else field_names[0],
                             size='mwu_qval' if 'mwu_qval' in field_names else field_names[0],
                             groups=group_names,
                             fields=field_names)
            de_result['data'] = de_result_data
            de_results.append(de_result)

            if n_genes > 0:
                field_use = None
                for field in ['mwu_qval', 'auroc', 't_qval']:
                    if field in field_names:
                        field_use = field
                        break
                if field_use is not None:
                    field_ascending = field_use != 'auroc'
                    for group_name in group_names:
                        fc_column = '{}:log2FC'.format(group_name)
                        name = '{}:{}'.format(group_name, field_name)

                        idx_up = de_result_df[fc_column] > 0
                        df_up = de_result_df.loc[idx_up].sort_values(by=[name, fc_column],
                                                                     ascending=[field_ascending, False])
                        features = df_up[:n_genes].index.values
                        marker_results.append(dict(category='markers', name=str(group_name), features=features))

        for key in dataset.obs_keys():
            if pd.api.types.is_categorical_dtype(dataset.obs[key]) or pd.api.types.is_bool_dtype(
                    dataset.obs[key]) or pd.api.types.is_object_dtype(dataset.obs[key]):
                obs_cat.append(key)
            else:
                obs.append(key)

        # if 'metagenes' in dataset.uns:
        #     metagenes = dataset.uns['metagenes']
        #     schema_dict['metafeatures'] = metagenes['var'].index

        for key in dataset.obs_keys():
            if pd.api.types.is_categorical_dtype(dataset.obs[key]) and len(dataset.obs[key]) < 1000:
                category_to_order[key] = dataset.obs[key].cat.categories

        # spatial_node = adata.uns['spatial'] if 'spatial' in adata.uns else None
        #
        # if spatial_node is not None:
        #     spatial_node_keys = list(spatial_node.keys())  # list of datasets
        #     if len(spatial_node_keys) == 1:
        #         spatial_node = spatial_node[spatial_node_keys[0]]

        images_node = dataset.uns.get('images',
                                      [])  # list of {type:image or meta_image, name:image name, image:path to image, spot_diameter:Number}
        image_names = list(map(lambda x: x['name'], images_node))

        for key in dataset.obsm_keys():
            dim = dataset.obsm[key].shape[1]
            if 1 < dim <= 3:
                embedding = dict(name=key, dimensions=dim)
                if dim == 2:
                    try:
                        image_index = image_names.index(key)
                        embedding['spatial'] = images_node[image_index]
                    except ValueError:
                        pass

                embeddings.append(embedding)
        meta_images = dataset.uns.get('meta_images', [])
        for meta_image in meta_images:
            embeddings.append(meta_image)

        for key in dataset.uns.keys():
            if key.endswith('_colors'):
                field = key[0:len(key) - len('_colors')]
                if field in dataset.obs:
                    colors = dataset.uns[key]
                    if pd.api.types.is_categorical_dtype(dataset.obs[field]):
                        categories = dataset.obs[field].cat.categories
                        if len(categories) == len(colors):
                            color_map = dict()
                            for i in range(len(categories)):
                                color_map[str(categories[i])] = colors[i]
                            field_to_value_to_color[field] = color_map
                    else:
                        print("Skipping colors for {}".format(key))
    schema_dict = {'version': '1.0.0'}
    schema_dict['results'] = de_results
    schema_dict['colors'] = field_to_value_to_color
    schema_dict['markers'] = marker_results
    schema_dict['embeddings'] = embeddings
    schema_dict['categoryOrder'] = category_to_order
    var_ids = []
    modules = None
    for dataset in datasets:
        if dataset.uns.get('data_type', '') == DataType.MODULE:
            modules = dataset.var if modules is None else pd.concat((modules, dataset.var))
        else:
            var_ids += dataset.var.index.to_list()
    if modules is not None:
        modules.index.name = 'id'
        schema_dict['modules'] = modules.reset_index().to_dict(orient='records')
    schema_dict['var'] = var_ids
    schema_dict['obs'] = obs
    schema_dict['obsCat'] = obs_cat

    schema_dict['shape'] = [datasets[0].shape[0], len(var_ids)]
    return schema_dict
