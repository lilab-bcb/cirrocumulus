import pandas as pd
import scipy.sparse

from cirrocumulus.dotplot_aggregator import DotPlotAggregator
from cirrocumulus.embedding_aggregator import EmbeddingAggregator, get_basis
from cirrocumulus.feature_aggregator import FeatureAggregator
from cirrocumulus.ids_aggregator import IdsAggregator
from cirrocumulus.unique_aggregator import UniqueAggregator


def apply_filter(df, data_filter):
    keep_expr = get_filter_expr(df, data_filter)
    return df[keep_expr] if keep_expr is not None else df


def get_filter_expr(df, data_filter):
    if df is None:
        raise ValueError('df is None')

    keep_expr = None
    if data_filter is not None:
        user_filters = data_filter.get('filters', [])
        combine_filters = data_filter.get('combine', 'and')
        for filter_obj in user_filters:
            field = filter_obj[0]
            op = filter_obj[1]
            value = filter_obj[2]
            if isinstance(field, dict):  # selection box

                selected_points_basis = get_basis(field['basis'], field.get('nbins'),
                    field.get('agg'), field.get('ndim', '2'), field.get('precomputed', False))

                if 'points' in value:
                    p = value['points']
                    field = selected_points_basis['full_name'] if selected_points_basis[
                                                                      'nbins'] is not None else 'index'
                    if field == 'index':
                        keep = df.index.isin(p)
                    else:
                        keep = df[field].isin(p)
                else:
                    keep = None
                    for p in value['path']:
                        if 'z' in p:  # 3d
                            selection_keep = \
                                (df[selected_points_basis['coordinate_columns'][0]] >= p['x']) & \
                                (df[selected_points_basis['coordinate_columns'][0]] <= p['x'] + p['width']) & \
                                (df[selected_points_basis['coordinate_columns'][1]] >= p['y']) & \
                                (df[selected_points_basis['coordinate_columns'][1]] <= p['y'] + p[
                                    'height']) & \
                                (df[selected_points_basis['coordinate_columns'][2]] >= p['z']) & \
                                (df[selected_points_basis['coordinate_columns'][2]] <= p['z'] + p['depth'])
                        else:
                            selection_keep = \
                                (df[selected_points_basis['coordinate_columns'][0]] >= p['x']) & \
                                (df[selected_points_basis['coordinate_columns'][0]] <= p['x'] + p['width']) & \
                                (df[selected_points_basis['coordinate_columns'][1]] >= p['y']) & \
                                (df[selected_points_basis['coordinate_columns'][1]] <= p['y'] + p['height'])

                    keep = selection_keep | keep if keep is not None else selection_keep

            else:
                series = df[field]
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


def split_measures(measures):
    obs_measures = []
    var_measures = []
    for i in range(len(measures)):
        name, key_type = get_var_name_type(measures[i])
        if key_type == 'X':
            var_measures.append(name)
        elif key_type == 'obs':
            obs_measures.append(name)
        else:
            raise ValueError('Unknown key type: ' + key_type)
    return obs_measures, var_measures


def handle_embedding(dataset_api, dataset, basis, measures=[], dimensions=[], quick=True):
    obs_measures, var_measures = split_measures(measures)
    if basis['precomputed']:
        result = precomputed_embedding(dataset_api, dataset, basis, obs_measures, var_measures, dimensions)
    else:
        count = '__count' in var_measures
        if count:
            var_measures.remove('__count')
        df = dataset_api.read(dataset=dataset, obs_keys=dimensions + obs_measures,
            var_keys=var_measures, basis=[basis])

        result = EmbeddingAggregator(obs_measures=obs_measures,
            var_measures=var_measures, dimensions=dimensions,
            count=count,
            nbins=basis['nbins'],
            basis=basis,
            agg_function=basis['agg'],
            quick=quick).execute(df)
    return result


def handle_grouped_stats(dataset_api, dataset, measures=[], dimensions=[]):
    # all dot plot measures are in X
    result = {}
    if dataset_api.has_precomputed_stats(dataset):
        result['dotplot'] = precomputed_grouped_stats(dataset_api, dataset, measures, dimensions)
    else:
        df = dataset_api.read(dataset=dataset, obs_keys=dimensions, var_keys=measures)
        result['dotplot'] = DotPlotAggregator(var_measures=measures,
            dimensions=dimensions).execute(df)
    return result


def handle_stats(dataset_api, dataset, measures=[], dimensions=[]):
    result = {}
    obs_measures, var_measures = split_measures(measures)
    if dataset_api.has_precomputed_stats(dataset):
        result['summary'] = precomputed_summary(dataset_api, dataset, obs_measures, var_measures,
            dimensions)
    else:
        df = dataset_api.read(dataset=dataset, obs_keys=obs_measures + dimensions, var_keys=var_measures)
        result['summary'] = FeatureAggregator(obs_measures, var_measures, dimensions).execute(df)
    return result


def handle_export_dataset_filters(dataset_api, dataset, data_filters):
    import json

    reformatted_filters = []
    filter_names = []
    for data_filter_obj in data_filters:
        filter_value = json.loads(data_filter_obj['value'])
        filter_names.append(data_filter_obj['name'])
        reformatted_filters.append(filter_value)

    df_info = get_data(dataset_api, dataset, measures=['obs/index'], data_filters=reformatted_filters)
    result_df = pd.DataFrame(index=df_info['df']['index'])
    for i in range(len(reformatted_filters)):
        data_filter = reformatted_filters[i]
        filter_name = filter_names[i]
        df_filtered = apply_filter(df_info['df'], data_filter)
        df = pd.DataFrame(index=df_filtered['index'])
        df[filter_name] = True
        result_df = result_df.join(df, rsuffix='r')
    result_df.fillna(False, inplace=True)
    return result_df.to_csv()


# def handle_diff_exp(dataset_api, dataset, data_filter, var_range):
#     df_info = get_data(dataset_api, dataset, [], [], [], [data_filter])
#     df = df_info['df']
#     mask = get_filter_expr(df, data_filter)
#     return dataset_api.diff_exp(dataset, mask, var_range)


def handle_selection(dataset_api, dataset, embeddings=[], measures=[], dimensions=[], data_filter=None, stats=True):
    selection_info = get_selected_data(dataset_api, dataset, embeddings=embeddings, measures=measures,
        dimensions=dimensions, data_filter=data_filter)
    basis_objs = selection_info['basis']
    obs_measures = selection_info['obs_measures']
    var_measures = selection_info['var_measures']
    df = selection_info['df']

    result = {}
    if len(basis_objs) > 0:
        result['coordinates'] = {}
        for basis_obj in basis_objs:
            result['coordinates'][basis_obj['full_name']] = UniqueAggregator(
                basis_obj['full_name'] if basis_obj['nbins'] is not None else 'index').execute(df)
    if stats:
        result['summary'] = FeatureAggregator(obs_measures, var_measures, dimensions).execute(df)
    if len(dimensions) > 0 and len(var_measures) > 0:
        result['dotplot'] = DotPlotAggregator(var_measures=var_measures,
            dimensions=dimensions).execute(df)
    result['count'] = df.shape[0]
    return result


def handle_selection_ids(dataset_api, dataset, data_filter):
    selection_info = get_selected_data(dataset_api, dataset, measures=['obs/index'], data_filter=data_filter)
    return {'ids': IdsAggregator().execute(selection_info['df'])}


def get_data(dataset_api, dataset, embeddings=[], measures=[], dimensions=[], data_filters=[]):
    obs_measures, var_measures = split_measures(measures)
    var_keys_filter = []
    obs_keys_filter = []
    selected_points_filter_basis_list = []
    for data_filter in data_filters:
        _var_keys_filter, _obs_keys_filter, selected_points_filter_basis_list = data_filter_keys(data_filter)
        selected_points_filter_basis_list += selected_points_filter_basis_list
        var_keys_filter += _var_keys_filter
        obs_keys_filter += _obs_keys_filter
    obs_keys = list(set(obs_measures + obs_keys_filter + dimensions))
    var_keys = list(set(var_measures + var_keys_filter))
    basis_objs = []
    basis_keys = set()
    for embedding in embeddings:
        basis_obj = get_basis(embedding['basis'], embedding.get('nbins'), embedding.get('agg'),
            embedding.get('ndim', '2'), embedding.get('precomputed', False))
        if basis_obj['full_name'] not in basis_keys:
            basis_objs.append(basis_obj)
            basis_keys.add(basis_obj['full_name'])
    for basis_obj in selected_points_filter_basis_list:
        if basis_obj['full_name'] not in basis_keys:
            basis_objs.append(basis_obj)
            basis_keys.add(basis_obj['full_name'])
    df = dataset_api.read(dataset=dataset, obs_keys=obs_keys, var_keys=var_keys, basis=basis_objs)
    for basis_obj in basis_objs:
        if not basis_obj['precomputed'] and basis_obj['nbins'] is not None:
            EmbeddingAggregator.convert_coords_to_bin(df=df,
                nbins=basis_obj['nbins'],
                bin_name=basis_obj['full_name'],
                coordinate_columns=basis_obj['coordinate_columns'])
    return {'df': df, 'basis': basis_objs, 'obs_measures': obs_measures, 'var_measures': var_measures,
            'dimensions': dimensions}


def get_selected_data(dataset_api, dataset, embeddings=[], measures=[], dimensions=[], data_filter=None):
    df_info = get_data(dataset_api, dataset, embeddings, measures, dimensions,
        [data_filter] if data_filter is not None else [])
    df = apply_filter(df_info['df'], data_filter)
    df_info['df'] = df
    return df_info


def data_filter_keys(data_filter):
    var_keys = set()
    obs_keys = set()
    basis_list = []
    if data_filter is not None:
        user_filters = data_filter.get('filters', [])

        for i in range(len(user_filters)):
            user_filter = user_filters[i]
            key = user_filter[0]
            if isinstance(key, dict):
                basis = get_basis(key['basis'], key.get('nbins'), key.get('agg'),
                    key.get('ndim', '2'), key.get('precomputed', False))
                basis_list.append(basis)
            else:
                name, key_type = get_var_name_type(key)
                user_filter[0] = name
                if key_type == 'X':
                    var_keys.add(name)
                else:
                    obs_keys.add(name)
    return list(var_keys), list(obs_keys), basis_list
