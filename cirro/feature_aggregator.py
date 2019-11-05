import pandas as pd
from natsort.natsort import natsorted


def sum(x):
    return x.values.sum(dtype='f8')


def fraction_expressed(x):
    return (x > 0).sum()


class FeatureAggregator:

    def __init__(self, measures, dimensions):
        self.measure_df = None
        self.dimension_value_counts = {}
        self.measures = measures
        self.dimensions = dimensions

    def __agg(self, df):
        cat_dict = {}
        grouped_df = None
        measure_df = None
        if (len(self.measures) + len(self.dimensions)) > 0:
            grouped_df = df.groupby('__selected')

        for column in self.dimensions:
            value_counts = grouped_df.agg({column: lambda x: x.value_counts()})
            value_counts = value_counts[value_counts > 0]
            cat_dict[column] = value_counts
        if len(self.measures) > 0:
            agg_func = {'__count': 'sum'}
            for column in self.measures:
                agg_func[column] = ['min', 'max', sum, fraction_expressed]
            measure_df = grouped_df.agg(agg_func)
        return measure_df, cat_dict

    def add(self, df):
        measure_df, dimension_value_counts = self.__agg(df)
        if measure_df is not None:
            self.measure_df = pd.concat(
                (self.measure_df,
                 measure_df)) if self.measure_df is not None else measure_df
            general_agg_dict = {'min': 'min', 'max': 'max', 'sum': 'sum', 'fraction_expressed': 'sum',
                                '__count': 'sum'}
            agg_dict = {}
            for column in self.measure_df:
                agg_dict[column] = general_agg_dict[column[1] if len(column) == 2 else column]
            self.measure_df = self.measure_df.groupby(
                self.measure_df.index).agg(agg_dict)
        for key in dimension_value_counts:
            existing_counts = self.dimension_value_counts.get(key, None)
            batch_counts = dimension_value_counts[key]
            existing_counts = pd.concat(
                (existing_counts,
                 batch_counts)) if existing_counts is not None else batch_counts
            self.dimension_value_counts[key] = existing_counts.groupby(existing_counts.index).agg('sum')


    def collect(self):
        # TODO add variance
        # divide sums by __count to get mean
        if self.measure_df is not None:
            df = self.measure_df

            fraction_expressed_column_names = []
            mean_column_names = []

            for name in df.columns.unique(0):
                mean_column_names.append((name, 'mean'))
                if name != '__count':
                    fraction_expressed_column_names.append((name, 'fraction_expressed'))

            count_values = df.iloc[:, df.columns.get_level_values(0) == '__count'].values
            mean = df.iloc[:, df.columns.get_level_values(1) == 'sum'].values / count_values
            fraction_expressed_array = df.iloc[:,
                                       df.columns.get_level_values(1) == 'fraction_expressed'].values / count_values
            for i in range(mean.shape[1]):
                df[mean_column_names[i]] = mean[:, i]
            for i in range(fraction_expressed_array.shape[1]):
                df[fraction_expressed_column_names[i]] = fraction_expressed_array[:, i]
            # remove count
            df = df.iloc[:, ~(df.columns.get_level_values(0) == '__count')]
            self.measure_df = df
        for key in self.dimension_value_counts:
            df = self.dimension_value_counts[key]
            sorted_categories = natsorted(df.index)
            df = df.loc[sorted_categories]
            self.dimension_value_counts[key] = df
        return self.measure_df, self.dimension_value_counts
