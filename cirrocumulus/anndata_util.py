import anndata
import numpy as np
import pandas as pd

DATA_TYPE_MODULE = 'module'
DATA_TYPE_UNS_KEY = 'data_type'
ADATA_MODULE_UNS_KEY = 'anndata_module'


def get_base(adata):
    base = None
    if 'log1p' in adata.uns and adata.uns['log1p']['base'] is not None:
        base = adata.uns['log1p'][base]
    return base


def adata_to_df(adata):
    df = pd.DataFrame(adata.X, index=adata.obs.index, columns=adata.var.index)
    for key in adata.layers.keys():
        df2 = pd.DataFrame(adata.layers[key], index=adata.obs.index.astype(str) + '-{}'.format(key),
                           columns=adata.var.index)
        df = pd.concat((df, df2), axis=0)

    df = df.T.join(adata.var)
    df.index.name = 'id'
    return df.reset_index()


def get_scanpy_marker_keys(dataset):
    marker_keys = []
    try:
        dataset.uns  # dataset can be AnnData or zarr group
    except AttributeError:
        return marker_keys
    from collections.abc import Mapping
    for key in dataset.uns.keys():
        rank_genes_groups = dataset.uns[key]
        if isinstance(rank_genes_groups, Mapping) and 'names' in rank_genes_groups and (
                'pvals' in rank_genes_groups or 'pvals_adj' in rank_genes_groups or 'scores' in rank_genes_groups) and len(
            rank_genes_groups['names'][0]) > 0 and not isinstance(rank_genes_groups['names'][0][0], bytes):
            marker_keys.append(key)
    return marker_keys


def get_pegasus_marker_keys(dataset):
    marker_keys = []
    try:
        dataset.varm  # dataset can be AnnData or zarr group
    except AttributeError:
        return marker_keys
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


def dataset_schema(dataset, n_features=10):
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

    prior_marker_results = dataset.uns.get('markers', [])
    if isinstance(prior_marker_results, str):
        import json
        prior_marker_results = json.loads(prior_marker_results)
    marker_results += prior_marker_results

    de_results = []  # array of dicts containing params logfoldchanges, pvals_adj, scores, names
    category_to_order = {}
    embeddings = []
    field_to_value_to_color = dict()  # field -> value -> color

    for scanpy_marker_key in get_scanpy_marker_keys(dataset):
        rank_genes_groups = dataset.uns[scanpy_marker_key]
        has_fc = 'logfoldchanges' in rank_genes_groups
        min_fold_change = 1
        params = rank_genes_groups['params']
        if not isinstance(params, dict):
            from anndata._io.zarr import read_attribute
            params = {k: read_attribute(params[k]) for k in params.keys()}

        # pts and pts_rest in later scanpy versions
        rank_genes_groups_keys = list(rank_genes_groups.keys())
        for k in ['params', 'names']:
            if k in rank_genes_groups_keys:
                rank_genes_groups_keys.remove(k)
        if 'pvals' in rank_genes_groups_keys and 'pvals_adj' in rank_genes_groups_keys:
            rank_genes_groups_keys.remove('pvals')
        category = '{} ({})'.format(params['groupby'], scanpy_marker_key)
        de_result_name = category
        de_result_df = None
        group_names = rank_genes_groups['names'].dtype.names
        de_result = dict(id='cirro-{}'.format(scanpy_marker_key),
                         type='de',
                         readonly=True,
                         groups=group_names,
                         fields=rank_genes_groups_keys,
                         name=de_result_name)
        for group_name in group_names:
            group_df = pd.DataFrame(index=rank_genes_groups['names'][group_name][...])
            group_df = group_df[group_df.index != 'nan']
            for rank_genes_groups_key in rank_genes_groups_keys:
                values = rank_genes_groups[rank_genes_groups_key][group_name][...]
                column_name = '{}:{}'.format(group_name, rank_genes_groups_key)
                group_df[column_name] = values

            if de_result_df is None:
                de_result_df = group_df
            else:
                de_result_df = de_result_df.join(group_df, how='outer')
            if n_features > 0:
                markers_df = group_df
                if has_fc:
                    markers_df = group_df[group_df['{}:logfoldchanges'.format(group_name)] > min_fold_change]
                if len(markers_df) > 0:
                    marker_results.append(
                        dict(category=de_result_name, name=str(group_name),
                             features=markers_df.index[:n_features]))

        if de_results_format == 'records':
            de_result_data = de_result_df.reset_index().to_dict(orient='records')
        else:
            de_result_data = dict(index=de_result_df.index[...])
            for c in de_result_df:
                de_result_data[c] = de_result_df[c]

        de_result['params'] = params
        de_result['data'] = de_result_data
        de_results.append(de_result)
    for pg_marker_key in get_pegasus_marker_keys(dataset):
        de_res = dataset.varm[pg_marker_key]
        key_name = pg_marker_key
        if pg_marker_key.startswith('rank_genes_'):
            key_name = pg_marker_key[len('rank_genes_'):]
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
        de_result = dict(id='cirro-{}'.format(pg_marker_key), type='de', name=key_name,
                         color='log2FC' if 'log2FC' in field_names else field_names[0],
                         size='mwu_qval' if 'mwu_qval' in field_names else field_names[0],
                         groups=group_names,
                         fields=field_names)
        de_result['data'] = de_result_data
        de_results.append(de_result)

        if n_features > 0:
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
                    features = df_up[:n_features].index.values

                    marker_results.append(dict(category=key_name, name=str(group_name), features=features))

    categories_node = dataset.obs['__categories'] if '__categories' in dataset.obs else None
    for key in dataset.obs.keys():
        if categories_node is not None and (key == '__categories' or key == 'index'):
            continue
        val = dataset.obs[key]
        if categories_node is not None and key in categories_node:
            categories = categories_node[key][...]
            ordered = categories_node[key].attrs.get('ordered', False)
            val = pd.Categorical.from_codes(val[...], categories, ordered=ordered)

        if pd.api.types.is_categorical_dtype(val) or pd.api.types.is_bool_dtype(
                val) or pd.api.types.is_object_dtype(val):
            obs_cat.append(key)
        else:
            obs.append(key)
        if pd.api.types.is_categorical_dtype(val) and len(val) < 1000:
            category_to_order[key] = dataset.obs[key].cat.categories
            categories = val.cat.categories
            color_field = key + '_colors'
            if color_field in dataset.uns:
                colors = dataset.uns[color_field][...]
                if len(categories) == len(colors):
                    color_map = dict()
                    for j in range(len(categories)):
                        color_map[str(categories[j])] = colors[j]
                    field_to_value_to_color[color_field] = color_map

    # spatial_node = adata.uns['spatial'] if 'spatial' in adata.uns else None
    #
    # if spatial_node is not None:
    #     spatial_node_keys = list(spatial_node.keys())  # list of datasets
    #     if len(spatial_node_keys) == 1:
    #         spatial_node = spatial_node[spatial_node_keys[0]]

    images_node = dataset.uns.get('images',
                                  [])  # list of {type:image or meta_image, name:image name, image:path to image, spot_diameter:Number}
    image_names = list(map(lambda x: x['name'], images_node))
    # layers = []
    # try:
    #     dataset.layers  # dataset can be AnnData or zarr group
    #     layers = list(dataset.layers.keys())
    #     # adata.list_keys()
    # except AttributeError:
    #     pass

    for key in dataset.obsm.keys():
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
    schema_dict = {'version': '1.0.0'}
    schema_dict['results'] = de_results
    schema_dict['colors'] = field_to_value_to_color
    schema_dict['markers'] = marker_results
    schema_dict['embeddings'] = embeddings
    schema_dict['categoryOrder'] = category_to_order

    var_df = dataset.var
    if not isinstance(var_df, pd.DataFrame):
        from anndata._io.zarr import read_attribute
        var_df = read_attribute(dataset.var)
    var_df.index.name = 'id'
    schema_dict['var'] = var_df.reset_index().to_dict(orient='records')
    modules_df = None
    if ADATA_MODULE_UNS_KEY in dataset.uns and isinstance(dataset.uns[ADATA_MODULE_UNS_KEY], anndata.AnnData):
        modules_df = dataset.uns[ADATA_MODULE_UNS_KEY].var
        #  if not isinstance(module_var, pd.DataFrame):
        #             from anndata._io.zarr import read_attribute
        #             module_var = read_attribute(module_var)
    if modules_df is not None:
        modules_df.index.name = 'id'
        schema_dict['modules'] = modules_df.reset_index().to_dict(orient='records')

    schema_dict['obs'] = obs
    schema_dict['obsCat'] = obs_cat
    shape = dataset.shape if isinstance(dataset, anndata.AnnData) else dataset.X.attrs.shape
    schema_dict['shape'] = shape
    return schema_dict
