import numpy as np


def mean_agg(x):
    return x.mean()


def sum_agg(x):
    return x.values.sum()


def max_agg(x):
    if hasattr(x, 'sparse'):
        if x.sparse.npoints == 0:
            return x.sparse.fill_value
        return np.max(x.values)
    return x.max()


def purity_agg(x):
    value_counts = x.value_counts(sort=False)
    largest = value_counts.nlargest(1)
    return largest[0] / value_counts.sum()


def mode_agg(x):
    return x.mode()[0]


def get_basis(basis, nbins=None, agg=None, dimensions=2, precomputed=False):
    if isinstance(dimensions, str):
        dimensions = int(dimensions)
    coordinate_columns = []
    for i in range(dimensions):
        coordinate_columns.append(basis + '_' + str(i + 1))
    full_name = basis + '_' + str(dimensions)
    if nbins is not None:
        full_name = full_name + '_' + str(nbins) + '_' + str(agg)
    return {'name': basis, 'dimensions': dimensions, 'coordinate_columns': coordinate_columns, 'nbins': nbins,
            'agg': agg, 'full_name': full_name, 'precomputed': precomputed}


class EmbeddingAggregator:

    def __init__(self, obs_measures, var_measures, dimensions, count, nbins, basis, agg_function='max', quick=False):
        self.nbins = nbins

        if agg_function is not None:
            if agg_function == 'max':
                agg_function = max_agg
            elif agg_function == 'mean':
                agg_function = mean_agg
            elif agg_function == 'sum':
                agg_function = sum_agg
            else:
                raise ValueError('Unknown agg')
        self.agg_function = agg_function
        self.obs_measures = obs_measures
        self.var_measures = var_measures
        self.dimensions = dimensions
        self.add_count = count
        self.basis = basis
        self.quick = quick  # do not compute purity.

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
        var_measures = self.var_measures
        obs_measures = self.obs_measures
        dimensions = self.dimensions
        add_count = self.add_count
        basis = self.basis
        nbins = self.nbins
        agg_function = self.agg_function
        compute_purity = not self.quick and nbins is not None
        if nbins is not None:

            EmbeddingAggregator.convert_coords_to_bin(df=df,
                nbins=nbins,
                coordinate_columns=basis['coordinate_columns'],
                bin_name=basis['full_name'])

            if add_count:
                df['__count'] = 1.0

            # bin level summary, coordinates have already been converted
            agg_dict = {}
            full_basis_name = basis['full_name']
            for column in basis['coordinate_columns']:
                agg_dict[column] = 'min'
            for column in obs_measures + var_measures:
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
            result['bins'] = agg_df.index

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
        for column in var_measures + obs_measures:
            series = agg_df[column]
            if compute_purity:
                series = series[series.columns[0]]
            is_sparse = hasattr(series, 'sparse')
            if is_sparse:
                series = series.sparse.to_dense()
            result['values'][column] = series
        for column in basis['coordinate_columns']:
            series = agg_df[column]
            if compute_purity:
                series = series[series.columns[0]]
            result['coordinates'][column] = series
        return result
