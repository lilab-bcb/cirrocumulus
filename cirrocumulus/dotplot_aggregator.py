import pandas as pd
from natsort import natsorted
from pandas import CategoricalDtype


class DotPlotAggregator:

    def __init__(self, var_measures, dimensions):
        self.var_measures = var_measures
        self.dimensions = dimensions

    def execute(self, df):
        results = []
        # {categories:[], name:'', values:[{name:'', fractionExpressed:0, mean:0}]}
        var_measures = self.var_measures
        dimensions = self.dimensions
        if len(var_measures) == 0 or len(dimensions) == 0:
            return results

        def mean(x):
            return x.mean()

        def fraction_expressed(x):
            return (x.values != 0).sum() / len(x)

        for dimension in dimensions:
            if pd.api.types.is_categorical_dtype(df[dimension]):
                if not df[dimension].dtype.ordered:
                    df[dimension] = df[dimension].astype(
                        CategoricalDtype(natsorted(df[dimension].dtype.categories), ordered=True))
                if len(df[dimension].dtype.categories) <= 1:
                    continue
            agg_result = df.groupby(dimension, observed=True).agg([mean, fraction_expressed])

            values = []
            dotplot_result = {'categories': agg_result.index, 'name': dimension, 'values': values}
            for name in var_measures:
                series = agg_result[name]
                is_sparse = hasattr(series, 'sparse')
                if is_sparse:
                    series = series.sparse.to_dense()
                values.append({'name': name,
                               'fractionExpressed': series['fraction_expressed'],
                               'mean': series['mean']})
            results.append(dotplot_result)
        return results
