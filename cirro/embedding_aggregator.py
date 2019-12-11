import numpy as np

from cirro.simple_data import SimpleData


def get_basis(basis):
    if basis is not None:
        embedding_ndim = 2
        if basis.endswith('3d'):
            basis = basis[0:len(basis) - 3]
            embedding_ndim = 3
        coordinate_columns = []
        for i in range(embedding_ndim):
            coordinate_columns.append(basis + '_' + str(i + 1))
        return {'name': basis, 'dimensions': embedding_ndim, 'coordinate_columns': coordinate_columns}


def agg_cat(x):
    value_counts = x.value_counts(sort=False)
    largest = value_counts.nlargest(1)
    purity = largest[0] / value_counts.sum()
    return largest.index[0], purity


class EmbeddingAggregator:

    def __init__(self, obs_measures, var_measures, dimensions, count, nbins, basis, agg_function):
        self.nbins = nbins
        self.agg_function = agg_function
        self.obs_measures = obs_measures
        self.var_measures = var_measures
        self.dimensions = dimensions
        self.count = count
        self.add_count = count or nbins is not None
        self.basis = basis

    @staticmethod
    def get_bin_level_agg_dict(measures, coordinate_columns, agg_function):
        agg_func = {}
        for column in measures:
            if column == '__count':
                agg_func[column] = 'sum'
            else:
                agg_func[column] = agg_function
        for column in coordinate_columns:
            agg_func[column] = 'min'
        return agg_func


    @staticmethod
    def convert_coords_to_bin(df, nbins, coordinate_columns, coordinate_column_to_range=None):
        # replace coordinates with bin, set df index to bin
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
            df.index = df[coordinate_columns[0]] * nbins + df[coordinate_columns[1]]
        else:
            df.index = df[coordinate_columns[2]] + nbins * (
                    df[coordinate_columns[1]] + nbins * df[coordinate_columns[0]])

    def execute(self, adata):
        result = {'coordinates': {}, 'values': {}}

        df = SimpleData.to_df(adata, self.obs_measures, self.var_measures, self.dimensions, self.basis)
        if self.add_count:
            df['__count'] = 1.0
        if self.nbins is not None:
            # bin level summary, coordinates have already been converted
            agg_dict = EmbeddingAggregator.get_bin_level_agg_dict(self.obs_measures + self.var_measures,
                self.basis['coordinate_columns'], self.agg_function)
            for column in self.dimensions:
                agg_dict[column] = agg_cat
            df = df.groupby(df.index).agg(agg_dict)
            result['bins'] = df.index.to_list()
        for column in self.basis['coordinate_columns']:
            result['coordinates'][column] = df[column].to_list()
        for column in self.var_measures + self.obs_measures:
            result['values'][column] = df[column].to_list()
        for column in self.dimensions:
            if self.nbins is not None:
                mode, purity = df[column].str
                result['values'][column] = dict(mode=mode.astype('category').to_list(), purity=purity.to_list())
            else:
                result['values'][column] = df[column].to_list()
        return result
