from natsort.natsort import natsorted

from cirrocumulus.simple_data import SimpleData


class FeatureAggregator:

    def __init__(self, obs_measures=[], var_measures=[], dimensions=[]):
        self.obs_measures = obs_measures
        self.var_measures = var_measures
        self.dimensions = dimensions

    @staticmethod
    def add_to_result(df, result):
        if df is not None:
            measure_df = df.T
            for feature in measure_df:
                values = measure_df[feature]
                feature_result = {'min': float(values.loc['min']),
                                  'max': float(values.loc['max']),
                                  'sum': float(values.loc['sum']),
                                  'mean': float(values.loc['mean'])}
                if 'numExpressed' in values.index:
                    feature_result['numExpressed'] = int(values.loc['numExpressed'])
                result[feature] = feature_result

    def execute(self, df):
        result = {}

        for column in self.dimensions:
            df_counts = df.agg({column: lambda x: x.value_counts(sort=False)})
            sorted_categories = natsorted(df_counts.index)
            df_counts = df_counts.iloc[df_counts.index.get_indexer_for(sorted_categories)]
            dimension_summary = {'categories': df_counts.index.to_list(),
                                 'counts': df_counts[column].to_list()}
            result[column] = dimension_summary

        if len(self.var_measures) > 0:
            FeatureAggregator.add_to_result(SimpleData.X_stats(df, self.var_measures), result)
        if len(self.obs_measures) > 0:
            FeatureAggregator.add_to_result(SimpleData.obs_stats(df, self.obs_measures), result)
        return result
