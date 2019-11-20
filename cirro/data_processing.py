import numpy as np
import pandas as pd

from cirro.dotplot_aggregator import DotPlotAggregator
from cirro.embedding_aggregator import EmbeddingAggregator
from cirro.feature_aggregator import FeatureAggregator
from cirro.ids_aggregator import IdsAggregator
from cirro.unique_index_aggregator import UniqueIndexAggregator


def filter_df(df, data_filter):
    if data_filter is not None:
        selected_points = data_filter.get('selected_points', None)

        keep_expr = None
        user_filters = data_filter.get('filters', [])
        for filter_obj in user_filters:
            field = filter_obj[0]
            op = filter_obj[1]
            value = filter_obj[2]
            if op == 'in':
                keep = (df[field].isin(value)).values
            elif op == '>':
                keep = (df[field] > value).values
            elif op == '=':
                keep = (df[field] == value).values
            elif op == '<':
                keep = (df[field] < value).values
            elif op == '!=':
                keep = (df[field] != value).values
            elif op == '>=':
                keep = (df[field] >= value).values
            elif op == '<=':
                keep = (df[field] <= value).values
            else:
                raise ValueError('Unknown filter')
            keep_expr = keep_expr & keep if keep_expr is not None else keep
        if selected_points is not None:
            keep = df.index.isin(selected_points)
            keep_expr = keep_expr & keep if keep_expr is not None else keep
        df = df[keep_expr]
    return df


def precomputed_summary(dataset_api, dataset, measures, dimensions):
    result = {}
    for dimension in dimensions:
        df = dataset_api.table(dataset, ['__index_level_0__', 'counts'], None, 'counts_' + dimension).to_pandas()
        result[dimension] = dict(categories=df.index.values.tolist(), counts=df['counts'].values.tolist())
    if len(measures) > 0:
        df = dataset_api.table(dataset, measures + ['__index_level_0__'], None, 'measure_summary').to_pandas()
        for column in measures:
            if column in df:
                values = df[column]
                result[column] = {'min': float(values.loc['min']), 'max': float(values.loc['max']),
                                  'sum': float(values.loc['sum']), 'mean': float(values.loc['mean']),
                                  'num_expressed': int(values.loc['num_expressed'])}
    return result


def precomputed_dotplot(dataset_api, dataset, measures, dimensions):
    result = []
    # columns in table are tuples of (feature, 'mean') or (feature, 'fraction_expressed')
    # {categories:[], name:'', values:[{name:'', fractionExpressed:0, mean:0}]}
    if len(measures) > 0 and len(dimensions) > 0:
        measure_columns = []
        for measure in measures:
            measure_columns.append(str((measure, 'mean')))
            measure_columns.append(str((measure, 'fraction_expressed')))
        for dimension in dimensions:
            # index name is dimension
            df = dataset_api.table(dataset, measure_columns + [dimension], None,
                                            'statistics_' + str(dimension)).to_pandas()

            values = []
            dimension_result = dict(categories=df.index.values.tolist(), name=dimension, values=values)
            result.append(dimension_result)
            for measure in measures:
                values.append(
                    dict(name=measure, fractionExpressed=df[(measure, 'fraction_expressed')].values.tolist(),
                        mean=df[(measure, 'mean')].values.tolist()))
    return result


def precomputed_embedding(dataset_api, dataset, basis, agg_function, nbins, measures,
                          dimensions):
    keys = dimensions + measures
    if len(keys) == 0:
        keys = ['__count']
    df = dataset_api.table(dataset, keys + ['__index_level_0__'], basis,
                                    basis['name'] + '_' + str(nbins) + '_' + str(agg_function)).to_pandas()
    result = {'coordinates': {}, 'values': {}, 'bins': df.index.values.tolist()}
    for column in basis['coordinate_columns']:
        result['coordinates'][column] = df[column].values.tolist()
    for column in df:
        if column not in basis['coordinate_columns']:
            result['values'][column] = df[column].values.tolist()
    return result


def process_data(dataset_api, dataset, return_types, basis=None, nbins=None, embedding_measures=[],
                 embedding_dimensions=[], dotplot_measures=[], dotplot_dimensions=[], summary_measures=[],
                 summary_dimensions=[], agg_function='max', data_filter=None):
    precomputed_dataset_summary = None
    if 'summary' in dataset:
        precomputed_dataset_summary = dataset['summary']  # {embeddings:[{nbins:500, agg:'max'}]}
    extra_keys = list()
    return_types = set(return_types)
    if 'ids' in return_types:
        extra_keys.append('index')
    if data_filter is not None:
        user_filters = data_filter.get('filters', [])
        selected_points = data_filter.get('selected_points', None)
        if selected_points is not None:
            data_filter['selected_points'] = np.array(selected_points)
        for filter_obj in user_filters:
            extra_keys.append(filter_obj[0])
    all_keys = list(set(
        extra_keys + embedding_dimensions + embedding_measures + dotplot_dimensions + dotplot_measures + summary_dimensions + summary_measures))
    result = {}
    precomputed_return_types = set(['embedding', 'summary', 'dotplot'])
    if precomputed_dataset_summary is not None and len(precomputed_return_types.intersection(return_types)) > 0:
        if 'embedding' in return_types:
            result['embedding'] = precomputed_embedding(dataset_api, dataset, basis, agg_function, nbins,
                embedding_measures, embedding_dimensions)
        if 'summary' in return_types:
            result['summary'] = precomputed_summary(dataset_api, dataset, summary_measures, summary_dimensions)
        if 'dotplot' in return_types:
            result['dotplot'] = precomputed_dotplot(dataset_api, dataset, dotplot_measures, dotplot_dimensions)
        return result

    coordinate_column_to_range = None
    add_bin = nbins is not None and (
            'ids' in return_types or 'selectionCoordinates' in return_types or 'embedding' in return_types or 'selectionSummary' in return_types)
    if nbins is not None:
        coordinate_column_to_range = dataset_api.statistics(dataset, keys=[], basis=basis)
    if 'summary' in return_types:
        result['summary'] = FeatureAggregator(summary_measures, summary_dimensions)
    if 'embedding' in return_types:  # uses all data
        result['embedding'] = EmbeddingAggregator(embedding_measures, embedding_dimensions,
            len(embedding_measures) + len(embedding_dimensions) == 0,
            nbins, basis, agg_function, None)
    if 'dotplot' in return_types:  # uses all data
        result['dotplot'] = DotPlotAggregator(dotplot_measures, dotplot_dimensions)
    if 'selectionSummary' in return_types:  # summary by selection and feature
        result['selectionSummary'] = FeatureAggregator(summary_measures, summary_dimensions)
    if 'selectionCoordinates' in return_types:  # return selected bins or indices
        result['selectionCoordinates'] = UniqueIndexAggregator()
    if 'ids' in return_types:  # uses selection only
        result['ids'] = IdsAggregator()

    nrows = 0
    for table in dataset_api.tables(dataset, all_keys, basis):  # full dataset
        df = table.to_pandas()
        df['__count'] = 1.0
        if 'ids' in return_types:
            df = df.reset_index()
        df.index = pd.RangeIndex(nrows, nrows + len(df))
        nrows += len(df)
        if add_bin:  # add bin before filtering since we might need to filter on bin
            EmbeddingAggregator.convert_coords_to_bin(df=df,
                nbins=nbins,
                coordinate_columns=basis['coordinate_columns'],
                coordinate_column_to_range=coordinate_column_to_range)
        df = filter_df(df, data_filter)
        for key in result:
            result[key].add(df)
    collected_results = {}
    for key in result:
        collected_results[key] = result[key].collect()
    return collected_results
