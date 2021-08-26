import numpy as np
import pandas as pd


def mean_agg(x):
    return x.mean()


def sparse_sum_agg(x):
    return x.sparse.fill_value if x.sparse.npoints == 0 else x.sum()


def sparse_max_agg(x):
    return x.sparse.fill_value if x.sparse.npoints == 0 else np.max(x.values.sp_values)


def purity_agg(x):
    value_counts = x.value_counts(sort=False)
    largest = value_counts.nlargest(1)
    return largest[0] / value_counts.sum()


def mode_agg(x):
    return x.mode()[0]


class EmbeddingAggregator:

    def __init__(self, measures, dimensions, nbins, basis, agg_function='max', quick=True,
                 coords=True):
        self.nbins = nbins
        self.agg_function = agg_function
        self.measures = measures
        self.dimensions = dimensions
        self.add_count = nbins is not None and len(measures) == 0 and len(dimensions) == 0
        self.basis = basis
        self.quick = quick  # do not compute purity.
        self.coords = coords  # return coordinates

    @staticmethod
    def convert_coords_to_bin(df, nbins, coordinate_columns, bin_name, coordinate_column_to_range=None):
        # replace coordinates with bin, set df[name] to bin
        for name in coordinate_columns:
            values = df[name].values
            if coordinate_column_to_range is not None:
                view_column_range = coordinate_column_to_range[name]
                column_min = view_column_range[0]
                column_max = view_column_range[1]
            else:
                column_min = values.min()
                column_max = values.max()
            df[name] = np.floor(np.interp(values, [column_min, column_max], [0, nbins - 1])).astype(int)
        if len(coordinate_columns) == 2:
            df[bin_name] = df[coordinate_columns[0]] * nbins + df[coordinate_columns[1]]
        else:
            df[bin_name] = df[coordinate_columns[2]] + nbins * (
                    df[coordinate_columns[1]] + nbins * df[coordinate_columns[0]])

    def execute(self, df):
        result = {'coordinates': {}, 'values': {}}
        measures = self.measures
        dimensions = self.dimensions
        add_count = self.add_count
        basis = self.basis
        if basis is not None:
            result['name'] = basis['name']
        nbins = self.nbins
        agg_function = self.agg_function
        compute_purity = not self.quick and nbins is not None
        if nbins is not None:
            if basis['full_name'] not in df:
                EmbeddingAggregator.convert_coords_to_bin(df=df, nbins=nbins,
                                                          coordinate_columns=basis['coordinate_columns'],
                                                          bin_name=basis['full_name'])

            if add_count:
                df['__count'] = 1.0

            # bin level summary, coordinates have already been converted
            agg_dict = {}
            full_basis_name = basis['full_name']
            for column in basis['coordinate_columns']:
                agg_dict[column] = 'min'
            for column in measures:
                if hasattr(df[column], 'sparse'):
                    if agg_function == 'max':
                        sparse_agg = sparse_max_agg
                    elif agg_function == 'mean':
                        sparse_agg = mean_agg
                    elif agg_function == 'sum':
                        sparse_agg = sparse_sum_agg
                    agg_dict[column] = sparse_agg
                else:
                    agg_dict[column] = agg_function

            for column in dimensions:
                agg_dict[column] = (mode_agg, purity_agg) if compute_purity else mode_agg
            if add_count:
                agg_dict['__count'] = 'sum'

            agg_df = df.groupby(full_basis_name).agg(agg_dict)

            if add_count:
                series = agg_df['__count']
                if compute_purity:
                    series = series[series.columns[0]]
                result['values']['__count'] = series
            if self.coords:
                result['coordinates']['bins'] = agg_df.index

        else:  # no binning
            if add_count:
                result['values']['__count'] = np.ones(len(df))
            agg_df = df
        for column in dimensions:
            series = agg_df[column]
            if nbins is not None:
                if compute_purity:
                    result['values'][column] = dict(value=series['mode_agg'], purity=series['purity_agg'])
                else:
                    result['values'][column] = dict(value=series)
            else:
                result['values'][column] = series
        for column in measures:
            series = agg_df[column]
            if compute_purity:
                series = series[series.columns[0]]
            is_sparse = hasattr(series, 'sparse')
            if is_sparse:
                # result['values'][column] = series.sparse.to_dense()
                result['values'][column] = dict(indices=series.values.sp_index.indices, values=series.values.sp_values)
            else:
                if pd.api.types.is_categorical_dtype(series):
                    result['values'][column] = dict(values=series.values, categories=series.cat.categories.values)
                else:
                    result['values'][column] = series

        if self.coords:
            for column in basis['coordinate_columns']:
                series = agg_df[column]
                if compute_purity:
                    series = series[series.columns[0]]
                result['coordinates'][column] = series
        return result
