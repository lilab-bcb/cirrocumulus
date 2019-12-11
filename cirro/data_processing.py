import numpy as np
import pandas as pd
import scipy.sparse

from cirro.dotplot_aggregator import DotPlotAggregator
from cirro.embedding_aggregator import EmbeddingAggregator
from cirro.feature_aggregator import FeatureAggregator
from cirro.ids_aggregator import IdsAggregator
from cirro.simple_data import SimpleData
from cirro.unique_index_aggregator import UniqueIndexAggregator


def filter_adata(adata, data_filter):
    if data_filter is not None:
        selected_points = data_filter.get('selected_points', None)
        keep_expr = None
        user_filters = data_filter.get('filters', [])
        for filter_obj in user_filters:
            field = filter_obj[0]
            op = filter_obj[1]
            value = filter_obj[2]
            if field in adata.obs:
                X = adata.obs[field]
            else:
                X = adata.X
                X = X[:, SimpleData.get_var_index(adata, field)] if len(X.shape) == 2 else X
            if op == 'in':
                keep = (X.isin(value))
            elif op == '>':
                keep = (X > value)
            elif op == '=':
                keep = (X == value)
            elif op == '<':
                keep = (X < value)
            elif op == '!=':
                keep = (X != value)
            elif op == '>=':
                keep = (X >= value)
            elif op == '<=':
                keep = (X <= value)
            else:
                raise ValueError('Unknown filter')
            if scipy.sparse.issparse(keep):
                keep = keep.toarray().flatten()
            if isinstance(keep, pd.Series):
                keep = keep.values
            keep_expr = keep_expr & keep if keep_expr is not None else keep

        if selected_points is not None:
            keep = adata.obs.index.isin(selected_points)
            keep_expr = keep_expr & keep if keep_expr is not None else keep

        adata = SimpleData.view(adata, keep_expr)
    return adata


def precomputed_summary(dataset_api, dataset, obs_measures, var_measures, dimensions):
    result = {}
    if '__count' in var_measures:
        var_measures.remove('__count')
    for dimension in dimensions:
        df = dataset_api.read_summarized(dataset, obs_keys=[dimension], path='counts')
        result[dimension] = dict(categories=df['index'].values.tolist(), counts=df['value'].values.tolist())

    if (len(obs_measures) + len(var_measures)) > 0:
        df = dataset_api.read_summarized(dataset, var_keys=var_measures, obs_keys=obs_measures, path='statistics',
            rename=True)
        for column in obs_measures + var_measures:
            result[column] = {'min': float(df['{}_min'.format(column)][0]),
                              'max': float(df['{}_max'.format(column)][0]),
                              'sum': float(df['{}_sum'.format(column)][0]),
                              'mean': float(df['{}_mean'.format(column)][0])}
            if '{}_num_expressed'.format(column) in df:
                result[column]['numExpressed'] = int(df['{}_num_expressed'.format(column)][0])
    return result


def precomputed_dotplot(dataset_api, dataset, obs_measures, var_measures, dimensions):
    result = []
    if (len(obs_measures) + len(var_measures)) > 0 and len(dimensions) > 0:
        for dimension in dimensions:
            df = dataset_api.read_summarized(dataset, index=True, rename=True, obs_keys=obs_measures,
                var_keys=var_measures, path='grouped_statistics/' + dimension)
            values = []
            dimension_result = dict(categories=df['index'].values.tolist(), name=dimension, values=values)
            result.append(dimension_result)
            for measure in obs_measures + var_measures:
                fraction_expressed = df[measure + '_fraction_expressed'].values.tolist()
                mean = df[measure + '_mean'].values.tolist()
                values.append(dict(name=measure, fractionExpressed=fraction_expressed, mean=mean))
    return result


def precomputed_embedding(dataset_api, dataset, basis, agg_function, nbins, obs_measures, var_measures,
                          dimensions):
    if basis is None or agg_function is None or nbins is None:
        raise ValueError('basis, agg, and nbins must be provided')

    if (len(obs_measures) + len(var_measures) + len(dimensions)) == 0:
        obs_measures = ['__count']

    path = 'obsm_summary/' + basis['name']
    df = dataset_api.read_summarized(dataset, obs_keys=obs_measures + dimensions, var_keys=var_measures,
        index=True, rename=True, path=path)
    result = {'coordinates': {}, 'values': {}, 'bins': df['index'].values.tolist()}

    for column in basis['coordinate_columns']:
        result['coordinates'][column] = df[column].values.tolist()

    for key in obs_measures + dimensions + var_measures:
        result['values'][key] = df[key + '_value'].values.tolist()
    return result


return_type_choices = set(['ids', 'selectionCoordinates', 'embedding', 'selectionSummary', 'summary', 'dotplot'])
precomputed_return_types = set(['embedding', 'summary', 'dotplot'])


def get_var_name_type(key):
    index = key.find('.')
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
            raise ValueError('Unknown key type')
    return obs_measures, var_measures


# measures start with X.name or obs.name to distinguish obs and X, dimensions must be obs
def process_data(dataset_api, dataset, return_types, basis=None, nbins=None, embedding_measures=[],
                 embedding_dimensions=[], dotplot_measures=[], dotplot_dimensions=[], summary_measures=[],
                 summary_dimensions=[], agg_function='max', data_filter=None):
    precomputed_dataset_summary = None
    if 'summary' in dataset:
        precomputed_dataset_summary = dataset['summary']  # {embeddings:[{nbins:500, agg:'max'}]}
        if basis is not None:
            basis['name'] = basis['name'] + '_' + str(nbins) + '_' + str(agg_function)
    return_types = set(return_types)
    obs_keys = set()
    var_keys = set()
    if 'ids' in return_types:
        obs_keys.add('id')
    embedding_obs_measures, embedding_var_measures = split_measures(embedding_measures)
    dotplot_obs_measures, dotplot_var_measures = split_measures(dotplot_measures)
    summary_obs_measures, summary_var_measures = split_measures(summary_measures)

    obs_keys.update(embedding_obs_measures)
    var_keys.update(embedding_var_measures)
    obs_keys.update(embedding_dimensions)

    obs_keys.update(dotplot_obs_measures)
    var_keys.update(dotplot_var_measures)
    obs_keys.update(dotplot_dimensions)

    obs_keys.update(summary_obs_measures)
    var_keys.update(summary_var_measures)
    obs_keys.update(summary_dimensions)

    var_keys.update(summary_var_measures)

    if data_filter is not None:
        user_filters = data_filter.get('filters', [])
        selected_points = data_filter.get('selected_points', None)
        if selected_points is not None:
            data_filter['selected_points'] = np.array(selected_points)
        for i in range(len(user_filters)):
            user_filter = user_filters[i]
            key = user_filter[0]
            name, key_type = get_var_name_type(key)
            user_filter[0] = name
            if key_type == 'X':
                var_keys.add(name)
            else:
                obs_keys.add(name)
    result = {}

    if precomputed_dataset_summary is not None and len(precomputed_return_types.intersection(return_types)) > 0:
        if 'embedding' in return_types:
            result['embedding'] = precomputed_embedding(dataset_api, dataset, basis, agg_function, nbins,
                embedding_obs_measures, embedding_var_measures, embedding_dimensions)
        if 'summary' in return_types:
            result['summary'] = precomputed_summary(dataset_api, dataset, summary_obs_measures, summary_var_measures,
                summary_dimensions)
        if 'dotplot' in return_types:
            result['dotplot'] = precomputed_dotplot(dataset_api, dataset, dotplot_obs_measures, dotplot_var_measures,
                dotplot_dimensions)
        return result

    add_bin = nbins is not None and precomputed_dataset_summary is None and (
            'ids' in return_types or 'selectionCoordinates' in return_types or 'embedding' in return_types or 'selectionSummary' in return_types)
    if 'summary' in return_types:
        result['summary'] = FeatureAggregator(summary_obs_measures, summary_var_measures, summary_dimensions)
    if 'embedding' in return_types:  # uses all data

        result['embedding'] = EmbeddingAggregator(obs_measures=embedding_obs_measures,
            var_measures=embedding_var_measures, dimensions=embedding_dimensions,
            count=len(embedding_measures) + len(embedding_dimensions) == 0,
            nbins=nbins, basis=basis, agg_function=agg_function)
    if 'dotplot' in return_types:  # uses all data
        result['dotplot'] = DotPlotAggregator(obs_measures=dotplot_obs_measures, var_measures=dotplot_var_measures,
            dimensions=dotplot_dimensions)
    if 'selectionSummary' in return_types:  # summary by selection and feature
        result['selectionSummary'] = FeatureAggregator(summary_obs_measures, summary_var_measures, summary_dimensions)
    if 'selectionCoordinates' in return_types:  # return selected bins or indices
        result['selectionCoordinates'] = UniqueIndexAggregator()
    if 'ids' in return_types:  # uses selection only
        result['ids'] = IdsAggregator()

    adata = dataset_api.read(dataset=dataset, obs_keys=list(obs_keys), var_keys=list(var_keys),
        basis=basis)  # TODO apply filters when reading

    adata.obs['__count'] = 1.0
    if add_bin:  # add bin before filtering since we might need to filter on bin
        EmbeddingAggregator.convert_coords_to_bin(df=adata.obs,
            nbins=nbins,
            coordinate_columns=basis['coordinate_columns'],
            coordinate_column_to_range=None)
    adata = filter_adata(adata, data_filter)

    for key in result:
        result[key] = result[key].execute(adata)
    return result
