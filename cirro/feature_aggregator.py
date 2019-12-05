import pandas as pd
from cirro.simple_data import SimpleData
from natsort.natsort import natsorted


class FeatureAggregator:

    def __init__(self, obs_measures, var_measures, dimensions):
        self.measure_df = None
        self.dimension_value_counts = {}
        self.obs_measures = obs_measures
        self.var_measures = var_measures
        self.dimensions = dimensions
        self.nrows = 0

    def __add(self, adata):
        dimension_dict = {}
        measure_df = None

        if adata.shape[0] > 0:
            for column in self.dimensions:
                value_counts = adata.obs.agg({column: lambda x: x.value_counts(sort=False)})
                value_counts = value_counts[value_counts[column] > 0]
                dimension_dict[column] = value_counts
            measure_df = None
            if len(self.var_measures) > 0:
                measure_df = SimpleData.X_stats(adata, self.var_measures)
            if len(self.obs_measures) > 0:
                obs_stats = SimpleData.obs_stats(adata, self.obs_measures)
                measure_df = obs_stats if measure_df is None else pd.concat((measure_df, obs_stats))

        return measure_df, dimension_dict

    def add(self, adata):
        self.nrows += adata.shape[0]
        measure_df, dimension_value_counts = self.__add(adata)
        if measure_df is not None:
            if self.measure_df is not None:
                self.measure_df = pd.concat((self.measure_df, measure_df))
                agg_dict = {'min': 'min', 'max': 'max', 'sum': 'sum', 'numExpressed': 'sum'}
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
                                   'numExpressed': int(values.loc['numExpressed'])}
        return result
