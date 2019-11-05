import numpy as np
import pandas as pd


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


class EmbeddingAggregator:

    def __init__(self, measures, dimensions, count, nbins, basis, agg_function):
        self.nbins = nbins
        self.agg_function = agg_function
        self.measures = measures
        self.dimensions = dimensions
        self.count = count
        if count or nbins is not None:
            self.measures = self.measures + ['__count']
        self.basis = basis
        self.measure_df = None
        self.binned_categorical_value_counts_dict = {}
        self.nrows = 0
        self.agg_dict = EmbeddingAggregator.get_bin_level_agg_dict(self.measures, self.basis['coordinate_columns'],
            'sum' if agg_function == 'mean' else agg_function)
        self.dimension_df = None


    @staticmethod
    def get_bin_level_agg_dict(measures, coordinate_columns, agg_function):
        agg_func = {}
        for column in measures:
            if column == '__count':
                agg_func[column] = 'sum'
            else:
                agg_func[column] = agg_function
        for column in coordinate_columns:
            agg_func[column] = lambda x: x.values[0]
        return agg_func


    @staticmethod
    def convert_coords_to_bin(df, nbins, coordinate_columns, coordinate_column_to_range):
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

    def collect(self):
        if self.nbins is not None and self.agg_function == 'mean':
            count = self.measure_df['__count']
            self.measure_df = self.measure_df.apply(
                lambda x: x / count if x.name != '__count' and x.name not in self.basis['coordinate_columns'] else x)
        if self.nbins is not None and len(self.dimensions) > 0:
            agg_dict = {}
            for column in self.basis['coordinate_columns']:
                agg_dict[column] = lambda x: x.values[0]
            for column in self.dimensions:
                agg_dict[column] = lambda x: x.mode()[0]
            self.dimension_df = self.dimension_df.groupby(self.dimension_df.index).agg(agg_dict)
        return self.measure_df, self.dimension_df

    def add(self, df):
        self.nrows += len(df)
        if self.nbins is not None:
            if len(self.dimensions) > 0:
                dimension_df = df[self.dimensions + self.basis['coordinate_columns']]
                self.dimension_df = pd.concat(
                    (self.dimension_df, dimension_df)) if self.dimension_df is not None else dimension_df
            # bin level summary, coordinates have already been converted
            df = df.groupby(df.index).agg(self.agg_dict)
        else:
            df = df[self.measures + self.dimensions + self.basis['coordinate_columns']]
        self.measure_df = pd.concat((self.measure_df, df)) if self.measure_df is not None else df
        if self.nbins is not None:
            self.measure_df = self.measure_df.groupby(self.measure_df.index).agg(self.agg_dict)
