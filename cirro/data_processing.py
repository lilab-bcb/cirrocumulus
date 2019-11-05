import numpy as np
import pandas as pd

from cirro.dotplot_aggregator import DotPlotAggregator
from cirro.embedding_aggregator import EmbeddingAggregator
from cirro.feature_aggregator import FeatureAggregator
from cirro.ids_aggregator import IdsAggregator
from cirro.selection_aggregator import SelectionAggregator


def mark_selected(df, data_filter):
    if data_filter is not None:
        selected_points = data_filter.get('selected_points', None) if data_filter is not None else None
        filters = []
        if selected_points is not None:
            filters.append(df.index.isin(selected_points))
        user_filters = data_filter.get('filters', [])
        for filter_obj in user_filters:
            field = filter_obj[0]
            op = filter_obj[1]
            value = filter_obj[2]
            if op == 'in':
                filters.append((df[field].isin(value)))
            elif op == '>':
                filters.append((df[field] > value))
            elif op == '=':
                filters.append((df[field] == value))
            elif op == '<':
                filters.append((df[field] < value))
            else:
                raise ValueError('Unknown filter')

        df['__selected'] = np.logical_and(*filters) if len(filters) > 1 else filters[0]
    else:
        df['__selected'] = True
    return df


def process_data(dataset_api, dataset, return_types, basis=None, nbins=None, embedding_measures=[],
                 embedding_dimensions=[], dotplot_measures=[], dotplot_dimensions=[], summary_measures=[],
                 summary_dimensions=[], agg_function='mean', data_filter=None):
    extra_keys = list()
    return_types = set(return_types)
    if 'ids' in return_types:
        extra_keys.append('index')
    if data_filter is not None:
        user_filters = data_filter.get('filters', [])
        for filter_obj in user_filters:
            extra_keys.append(filter_obj[0])
    all_keys = list(set(
        extra_keys + embedding_dimensions + embedding_measures + dotplot_dimensions + dotplot_measures + summary_dimensions + summary_measures))
    result = {}
    coordinate_column_to_range = None

    if nbins is not None:
        coordinate_column_to_range = dataset_api.statistics(dataset, basis['coordinate_columns'])
    if 'summary' in return_types:  # summary by selection and feature
        result['summary'] = FeatureAggregator(summary_measures, summary_dimensions)
    if 'embedding' in return_types:  # uses all data
        result['embedding'] = EmbeddingAggregator(embedding_measures, embedding_dimensions,
            len(embedding_measures) + len(embedding_dimensions) == 0,
            nbins, basis, agg_function)
    if 'dotplot' in return_types:  # uses all data
        result['dotplot'] = DotPlotAggregator(dotplot_measures, dotplot_dimensions)
    if 'selection' in return_types:  # return selected bins or indices
        result['selection'] = SelectionAggregator()
    if 'ids' in return_types:  # uses selection only
        result['ids'] = IdsAggregator()

    nrows = 0
    add_bin = nbins is not None and (
            'ids' in return_types or 'selection' in return_types or 'embedding' in return_types)
    for table in dataset_api.tables(dataset, all_keys, basis):
        df = table.to_pandas()
        df['__count'] = 1.0
        if 'ids' in return_types:
            df = df.reset_index()
        df.index = pd.RangeIndex(nrows, nrows + len(df))
        nrows += len(df)
        if add_bin:
            EmbeddingAggregator.convert_coords_to_bin(df=df,
                nbins=nbins,
                coordinate_columns=basis['coordinate_columns'],
                coordinate_column_to_range=coordinate_column_to_range)
        mark_selected(df, data_filter)
        for key in result:
            result[key].add(df)
    return result
