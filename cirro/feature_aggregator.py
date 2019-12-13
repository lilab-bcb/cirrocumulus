import pandas as pd
from natsort.natsort import natsorted

from cirro.simple_data import SimpleData


class FeatureAggregator:

    def __init__(self, obs_measures, var_measures, dimensions):
        self.obs_measures = obs_measures
        self.var_measures = var_measures
        self.dimensions = dimensions

    def execute(self, adata):
        result = {}
        if adata.shape[0] > 0:
            for column in self.dimensions:
                df = adata.obs.agg({column: lambda x: x.value_counts(sort=False)})
                sorted_categories = natsorted(df.index)
                df = df.loc[sorted_categories]
                dimension_summary = {'categories': df.index.to_list(),
                                     'counts': df[column].to_list()}
                result[column] = dimension_summary
            measure_df = None
            if len(self.var_measures) > 0:
                measure_df = SimpleData.X_stats(adata, self.var_measures)
            if len(self.obs_measures) > 0:
                obs_stats = SimpleData.obs_stats(adata, self.obs_measures)
                measure_df = obs_stats if measure_df is None else pd.concat((measure_df, obs_stats))
            if measure_df is not None:
                measure_df = measure_df.T
                for feature in measure_df:
                    values = measure_df[feature]
                    feature_result = {'min': float(values.loc['min']),
                                      'max': float(values.loc['max']),
                                      'sum': float(values.loc['sum']),
                                      'mean': float(values.loc['mean'])}
                    result[feature] = feature_result
                    if 'numExpressed' in values.index:
                        feature_result['numExpressed'] = int(values.loc['numExpressed'])
        return result
