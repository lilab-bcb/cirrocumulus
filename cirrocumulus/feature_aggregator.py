from cirrocumulus.anndata_util import X_stats, obs_stats


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

    def execute(self, adata):
        result = {}
        for column in self.dimensions:
            df_counts = adata.obs.agg({column: lambda x: x.value_counts(sort=False)})
            dimension_summary = {'categories': df_counts.index,
                                 'counts': df_counts[column]}
            result[column] = dimension_summary
        if len(self.var_measures) > 0:
            FeatureAggregator.add_to_result(X_stats(adata), result)
        if len(self.obs_measures) > 0:
            FeatureAggregator.add_to_result(obs_stats(adata, self.obs_measures), result)
        return result
