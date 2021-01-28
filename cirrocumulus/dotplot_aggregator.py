import pandas as pd


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

        for d in dimensions:
            dimension_name = d
            if isinstance(d, list):
                if len(d) > 1:
                    dimension_name = '-'.join(d)
                    df[dimension_name] = df[d[0]].astype(str).str.cat(df[d[1:]].astype(str), sep="-").astype('category')
                else:
                    dimension_name = d[0]
            if pd.api.types.is_categorical_dtype(df[dimension_name]) and len(df[dimension_name].dtype.categories) <= 1:
                continue
            agg_result = df.groupby(dimension_name, observed=True).agg([mean, fraction_expressed])

            values = []
            dotplot_result = {'categories': agg_result.index, 'name': dimension_name, 'values': values}
            for var_measure in var_measures:
                series = agg_result[var_measure]
                is_sparse = hasattr(series, 'sparse')
                if is_sparse:
                    series = series.sparse.to_dense()
                values.append({'name': var_measure,
                               'percentExpressed': 100 * series['fraction_expressed'],
                               'mean': series['mean']})
            results.append(dotplot_result)
        return results
