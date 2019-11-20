import pandas as pd
from natsort.natsort import natsorted


def sum(x):
    return x.values.sum(dtype='f8')


def num_expressed(x):
    return (x > 0).sum()


class FeatureAggregator:

    def __init__(self, measures, dimensions):
        self.measure_df = None
        self.dimension_value_counts = {}
        self.measures = measures
        self.dimensions = dimensions
        self.nrows = 0

    def __add(self, df):
        dimension_dict = {}
        measure_df = None

        if len(df) > 0:
            for column in self.dimensions:
                value_counts = df.agg({column: lambda x: x.value_counts(sort=False)})
                value_counts = value_counts[value_counts[column] > 0]
                dimension_dict[column] = value_counts
            if len(self.measures) > 0:
                measure_df = df[self.measures].agg(['min', 'max', sum, num_expressed]).T
                # measures on columns, stats on rows, transpose so that stats are on columns
        return measure_df, dimension_dict

    def add(self, df):
        self.nrows += len(df)
        df = df[self.dimensions + self.measures]
        measure_df, dimension_value_counts = self.__add(df)
        if measure_df is not None:
            if self.measure_df is not None:
                self.measure_df = pd.concat((self.measure_df, measure_df))
                agg_dict = {'min': 'min', 'max': 'max', 'sum': 'sum', 'num_expressed': 'sum'}
                self.measure_df = self.measure_df.groupby(self.measure_df.index).agg(agg_dict)
            else:
                self.measure_df = measure_df
        for key in dimension_value_counts:
            existing_counts = self.dimension_value_counts.get(key, None)
            batch_counts = dimension_value_counts[key]
            if existing_counts is not None:
                existing_counts = pd.concat((existing_counts, batch_counts))
                self.dimension_value_counts[key] = existing_counts.groupby(existing_counts.index).agg('sum')
            else:
                self.dimension_value_counts[key] = batch_counts

    def collect(self):
        # TODO add variance
        # we return sum so that we can compute unselected mean using count and sum
        result = {}
        for key in self.dimension_value_counts:
            df = self.dimension_value_counts[key]
            sorted_categories = natsorted(df.index)
            df = df.loc[sorted_categories]
            dimension_summary = {'categories': df.index.values.tolist(),
                                 'counts': df[key].values.tolist()}
            result[key] = dimension_summary

        if self.measure_df is not None:
            df = self.measure_df.T
            for feature in df:
                values = df[feature]
                result[feature] = {'min': float(values.loc['min']),
                                   'max': float(values.loc['max']),
                                   'sum': float(values.loc['sum']),
                                   'mean': float(values.loc['sum'] / self.nrows),
                                   'num_expressed': int(values.loc['num_expressed'])}
        return result
