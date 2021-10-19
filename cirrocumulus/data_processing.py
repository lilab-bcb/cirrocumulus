import pandas as pd
import scipy.sparse

from cirrocumulus.anndata_util import ADATA_MODULE_UNS_KEY
from cirrocumulus.dotplot_aggregator import DotPlotAggregator
from cirrocumulus.feature_aggregator import FeatureAggregator
from cirrocumulus.ids_aggregator import IdsAggregator
from cirrocumulus.unique_aggregator import UniqueAggregator


def get_mask(dataset_api, dataset, data_filters):
    measures = set()
    dimensions = set()
    basis = set()
    for data_filter in data_filters:
        _measures, _dimensions, _basis = data_filter_keys(data_filter)
        measures.update(_measures)
        dimensions.update(_dimensions)
        basis.update(_basis)
    keys = get_type_to_measures(measures)
    keys['obs'] = list(dimensions)
    keys['basis'] = list(basis)
    adata = dataset_api.read_dataset(keys=keys, dataset=dataset)
    result = []
    for data_filter in data_filters:
        result.append(get_filter_expr(adata, data_filter))
    return result, adata


def apply_filter(adata, data_filter):
    keep_expr = get_filter_expr(adata, data_filter)
    return adata[keep_expr] if keep_expr is not None else adata


def get_filter_str(data_filter):
    user_filters = data_filter.get('filters', [])
    combine_filters = data_filter.get('combine', 'and')
    s = []
    for filter_obj in user_filters:
        field = filter_obj[0]
        op = filter_obj[1]
        value = filter_obj[2]
        if op == 'in':
            op = ''
        if not isinstance(field, dict) and not field == '__index':
            if isinstance(value, list):
                value = '_'.join(value)
                s.append(field + '-' + value)
            else:
                return None
        else:
            return None
    return ''.join(s).replace('/', '-').replace(' ', '-').replace('\\', '-')


def get_filter_expr(adata, data_filter):
    keep_expr = None
    if data_filter is not None:
        user_filters = data_filter.get('filters', [])
        combine_filters = data_filter.get('combine', 'and')
        for filter_obj in user_filters:
            field = filter_obj[0]
            op = filter_obj[1]
            value = filter_obj[2]
            if isinstance(field, dict):  # selection box
                # selected_points_basis = get_basis(field['basis'], field.get('nbins'),
                #                                   field.get('agg'), field.get('ndim', '2'),
                #                                   field.get('precomputed', False))

                if 'points' in value:
                    p = value['points']
                    # list of indices
                    keep = pd.RangeIndex(adata.shape[0]).isin(p)

                # else:
                #     keep = None
                #     for p in value['path']:
                #         if 'z' in p:  # 3d
                #             selection_keep = \
                #                 (df[selected_points_basis['coordinate_columns'][0]] >= p['x']) & \
                #                 (df[selected_points_basis['coordinate_columns'][0]] <= p['x'] + p['width']) & \
                #                 (df[selected_points_basis['coordinate_columns'][1]] >= p['y']) & \
                #                 (df[selected_points_basis['coordinate_columns'][1]] <= p['y'] + p[
                #                     'height']) & \
                #                 (df[selected_points_basis['coordinate_columns'][2]] >= p['z']) & \
                #                 (df[selected_points_basis['coordinate_columns'][2]] <= p['z'] + p['depth'])
                #         else:
                #             selection_keep = \
                #                 (df[selected_points_basis['coordinate_columns'][0]] >= p['x']) & \
                #                 (df[selected_points_basis['coordinate_columns'][0]] <= p['x'] + p['width']) & \
                #                 (df[selected_points_basis['coordinate_columns'][1]] >= p['y']) & \
                #                 (df[selected_points_basis['coordinate_columns'][1]] <= p['y'] + p['height'])

                # keep = selection_keep | keep if keep is not None else selection_keep

            elif field == '__index':
                import numpy as np
                keep = np.zeros(adata.shape[0], dtype=bool)
                keep[value] = True
            else:
                series = adata.obs[field] if field in adata.obs else adata[:, adata.var.index == field].X[:, 0]
                if op == 'in':
                    keep = (series.isin(value)).values
                elif op == '>':
                    keep = (series > value)
                elif op == '=':
                    keep = (series == value)
                elif op == '<':
                    keep = (series < value)
                elif op == '!=':
                    keep = (series != value)
                elif op == '>=':
                    keep = (series >= value)
                elif op == '<=':
                    keep = (series <= value)
                else:
                    raise ValueError('Unknown filter')

            if scipy.sparse.issparse(keep):
                keep = keep.toarray().flatten()
            if hasattr(keep, 'sparse'):
                keep = keep.sparse.to_dense()
            if isinstance(keep, pd.Series):
                keep = keep.values
            if keep_expr is not None:
                if combine_filters == 'and':
                    keep_expr = keep_expr & keep
                else:
                    keep_expr = keep_expr | keep
            else:
                keep_expr = keep

    return keep_expr


def precomputed_summary(dataset_api, dataset, obs_measures, var_measures, dimensions):
    if '__count' in var_measures:
        var_measures.remove('__count')
    return dataset_api.read_precomputed_stats(dataset, obs_keys=dimensions + obs_measures, var_keys=var_measures)


def precomputed_grouped_stats(dataset_api, dataset, var_measures, dimensions):
    if (len(var_measures)) > 0 and len(dimensions) > 0:
        return dataset_api.read_precomputed_grouped_stats(dataset,
                                                          var_keys=var_measures, obs_keys=dimensions)
    return []


def precomputed_embedding(dataset_api, dataset, basis, obs_measures, var_measures,
                          dimensions):
    if (len(obs_measures) + len(var_measures) + len(dimensions)) == 0:
        obs_measures = ['__count']
    return dataset_api.read_precomputed_basis(dataset, obs_keys=obs_measures + dimensions, var_keys=var_measures,
                                              basis=basis)


def get_var_name_type(key):
    index = key.find('/')
    if index == -1:
        return key, 'X'
    else:
        key_type = key[0:index]
        name = key[index + 1:]
        return name, key_type


def get_type_to_measures(measures):
    type2measures = dict(X=[], obs=[])
    for measure in measures:
        name, key_type = get_var_name_type(measure)
        type_measures = type2measures.get(key_type)
        if type_measures is None:
            type_measures = []
            type2measures[key_type] = type_measures
        type_measures.append(name)
    return type2measures


def check_bin_input(nbins):
    if nbins is not None:
        nbins = int(nbins)
        nbins = min(1000, nbins)
        if nbins <= 0:
            nbins = None
    return nbins


def handle_export_dataset_filters(dataset_api, dataset, data_filters):
    import json

    reformatted_filters = []
    filter_names = []
    for data_filter_obj in data_filters:
        filter_value = json.loads(data_filter_obj['value'])
        filter_names.append(data_filter_obj['name'])
        reformatted_filters.append(filter_value)

    df = get_adata(dataset_api, dataset, measures=['obs/index'], data_filters=reformatted_filters)
    result_df = pd.DataFrame(index=df['index'])
    for i in range(len(reformatted_filters)):
        data_filter = reformatted_filters[i]
        filter_name = filter_names[i]
        df_filtered = apply_filter(df, data_filter)
        df = pd.DataFrame(index=df_filtered['index'])
        df[filter_name] = True
        result_df = result_df.join(df, rsuffix='r')
    result_df.fillna(False, inplace=True)
    return result_df.to_csv()


# embedding - list of basis and coords (whether to return coordinates). For binned embeddings, also measures and dimensions.
def handle_data(dataset_api, dataset, embedding_list=None, values=None, grouped_stats=None, stats=None, selection=None):
    dimensions = set()
    measures = set()
    basis_keys = set()
    if embedding_list is not None:
        for embedding in embedding_list:
            basis_keys.add(embedding['name'])

    if values is not None:
        dimensions.update(values.get('dimensions', []))
        measures.update(values.get('measures', []))

    if selection is not None:
        data_filter = selection.get('filter')
        var_keys_filter, obs_keys_filter, selected_points_filter_basis_list = data_filter_keys(data_filter)
        selection['basis'] = selected_points_filter_basis_list
        measures.update(var_keys_filter)
        dimensions.update(obs_keys_filter)
        dimensions.update(selection.get('dimensions', []))
        measures.update(selection.get('measures', []))
        selection_embeddings = selection.get('embeddings', [])

        for embedding in selected_points_filter_basis_list + selection_embeddings:
            basis_keys.add(embedding['name'])

    if grouped_stats is not None:
        grouped_stats_dimensions = grouped_stats.get('dimensions', [])
        for d in grouped_stats_dimensions:
            if isinstance(d, list):
                dimensions.update(d)
            else:
                dimensions.add(d)
        measures.update(grouped_stats.get('measures', []))
    if stats is not None:
        dimensions.update(stats.get('dimensions', []))
        measures.update(stats.get('measures', []))

    keys = get_type_to_measures(measures)
    keys['obs'] += list(dimensions)
    keys['basis'] = list(basis_keys)
    adata = dataset_api.read_dataset(dataset=dataset, keys=keys)
    results = {}
    if values is not None:
        dimensions = values.get('dimensions', [])
        measures = values.get('measures', [])
        type2measures = get_type_to_measures(measures)
        results['values'] = {}
        for key in type2measures['obs'] + dimensions:
            series = adata.obs[key]
            results['values'][key] = series
            if pd.api.types.is_categorical_dtype(series):
                results['values'][key] = dict(values=series.values.codes, categories=series.cat.categories.values)
            else:
                results['values'][key] = series

        if adata.uns.get(ADATA_MODULE_UNS_KEY) is not None:
            adata_modules = adata.uns[ADATA_MODULE_UNS_KEY]
            for i in range(len(adata_modules.var.index)):
                x = adata_modules.X[:, i]
                results['values'][adata_modules.var.index[i]] = x
        if adata.X is not None:
            is_sparse = scipy.sparse.issparse(adata.X)
            for i in range(len(adata.var.index)):
                x = adata.X[:, i]
                if is_sparse:
                    indices = x.indices
                    data = x.data
                    results['values'][adata.var.index[i]] = dict(indices=indices, values=data)
                else:
                    results['values'][adata.var.index[i]] = x

    if embedding_list is not None and len(embedding_list) > 0:
        results['embeddings'] = []
        for key in adata.obsm.keys():
            m = adata.obsm[key]
            ndim = m.shape[1]
            coordinates = dict()
            embedding = dict(name=key, dimensions=ndim, coordinates=coordinates)
            results['embeddings'].append(embedding)
            for i in range(ndim):
                coordinates['{}_{}'.format(key, i + 1)] = m[:, i]

    if grouped_stats is not None:
        results['distribution'] = DotPlotAggregator(var_measures=grouped_stats.get('measures', []),
                                                    dimensions=grouped_stats.get('dimensions', [])).execute(adata)
    if stats is not None:
        dimensions = stats.get('dimensions', [])
        measures = stats.get('measures', [])
        type2measures = get_type_to_measures(measures)
        results['summary'] = FeatureAggregator(type2measures['obs'], type2measures['X'],
                                               dimensions).execute(adata)
    if selection is not None:
        results['selection'] = {}
        dimensions = selection.get('dimensions', [])
        measures = selection.get('measures', [])
        type2measures = get_type_to_measures(measures)
        # basis_list = selection.get('basis', [])
        selection_embeddings = selection.get('embeddings', [])
        df = apply_filter(adata, data_filter)
        if len(selection_embeddings) > 0:
            results['selection']['coordinates'] = {}
            for embedding in selection_embeddings:
                results['selection']['coordinates'][embedding['name']] = UniqueAggregator('index').execute(df)
        var_measures = type2measures['X']
        results['selection']['summary'] = FeatureAggregator(type2measures['obs'], type2measures['X'],
                                                            dimensions).execute(df)
        if len(dimensions) > 0 and len(var_measures) > 0:
            results['selection']['distribution'] = DotPlotAggregator(var_measures=var_measures,
                                                                     dimensions=[dimensions]).execute(df)
        results['selection']['count'] = df.shape[0]

    return results


def handle_selection_ids(dataset_api, dataset, data_filter):
    df = get_selected_data(dataset_api, dataset, measures=['obs/index'], data_filter=data_filter)
    return {'ids': IdsAggregator().execute(df)}


def get_adata(dataset_api, dataset, embeddings=[], measures=[], dimensions=[], data_filters=[]):
    type2measures = get_type_to_measures(measures)
    var_keys_filter = []
    obs_keys_filter = []
    selected_points_filter_basis_list = []
    for data_filter in data_filters:
        _var_keys_filter, _obs_keys_filter, selected_points_filter_basis_list = data_filter_keys(data_filter)
        selected_points_filter_basis_list += selected_points_filter_basis_list
        var_keys_filter += _var_keys_filter
        obs_keys_filter += _obs_keys_filter
    obs_keys = list(set(type2measures['obs'] + obs_keys_filter + dimensions))
    var_keys = list(set(type2measures['X'] + var_keys_filter))
    basis_keys = set()
    for embedding in embeddings:
        basis_keys.add(embedding['name'])
    for embedding in selected_points_filter_basis_list:
        basis_keys.add(embedding['name'])

    adata = dataset_api.read_dataset(dataset=dataset, keys=dict(obs=obs_keys, X=var_keys, basis=list(basis_keys)))
    # for basis_obj in basis_objs:
    #     if not basis_obj['precomputed'] and basis_obj['nbins'] is not None:
    #         EmbeddingAggregator.convert_coords_to_bin(df=df,
    #                                                   nbins=basis_obj['nbins'],
    #                                                   bin_name=basis_obj['full_name'],
    #                                                   coordinate_columns=basis_obj['coordinate_columns'])
    return adata


def get_selected_data(dataset_api, dataset, embeddings=[], measures=[], dimensions=[], data_filter=None):
    df = get_adata(dataset_api, dataset, embeddings, measures, dimensions,
                   [data_filter] if data_filter is not None else [])
    df = apply_filter(df, data_filter)
    return df


def data_filter_keys(data_filter):
    var_keys = set()
    obs_keys = set()
    basis_keys = set()
    if data_filter is not None:
        user_filters = data_filter.get('filters', [])

        for i in range(len(user_filters)):
            user_filter = user_filters[i]
            key = user_filter[0]
            if isinstance(key, dict):
                basis_keys.add(key['name'])
            elif key == '__index':
                continue
            else:
                name, key_type = get_var_name_type(key)
                user_filter[0] = name
                if key_type == 'X':
                    var_keys.add(name)
                else:
                    obs_keys.add(name)
    return list(var_keys), list(obs_keys), list(basis_keys)
