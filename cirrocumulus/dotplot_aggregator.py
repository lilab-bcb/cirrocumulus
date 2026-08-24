import numpy as np
import pandas as pd
import scipy.sparse
from pandas import CategoricalDtype


class DotPlotAggregator:
    def __init__(self, var_measures, dimensions):
        self.var_measures = var_measures
        self.dimensions = dimensions

    def _to_frame(self, adata):
        """Collect the requested dimensions and measures into one DataFrame."""
        names = []
        for d in self.dimensions:
            names.extend(d if isinstance(d, list) else [d])

        columns = {name: adata.obs[name] for name in dict.fromkeys(names) if name in adata.obs}

        var_index = list(adata.var.index)
        for measure in self.var_measures:
            name = measure.split("/", 1)[-1]
            if name in var_index:
                values = adata.X[:, var_index.index(name)]
                if scipy.sparse.issparse(values):
                    values = values.todense()
                columns[measure] = np.asarray(values).ravel()
            elif name in adata.obs:
                columns[measure] = adata.obs[name]
        return pd.DataFrame(columns, index=adata.obs.index)

    def execute(self, adata):
        results = []
        # {categories:[], name:'', values:[{name:'', percentExpressed:0, mean:0}]}
        var_measures = self.var_measures
        dimensions = self.dimensions
        if len(var_measures) == 0 or len(dimensions) == 0:
            return results

        df = self._to_frame(adata)

        def mean(x):
            return x.mean()

        def percent_expressed(x):
            return 100 * ((x.values != 0).sum() / len(x))

        for d in dimensions:
            dimension_name = d
            if isinstance(d, list):
                if len(d) > 1:
                    dimension_name = "-".join(d)
                    df[dimension_name] = (
                        df[d[0]]
                        .astype(str)
                        .str.cat(df[d[1:]].astype(str), sep="-")
                        .astype("category")
                    )
                else:
                    dimension_name = d[0]
            if (
                isinstance(df[dimension_name].dtype, CategoricalDtype)
                and len(df[dimension_name].dtype.categories) <= 1
            ):
                continue
            # aggregate only the measures: with several dimensions the frame also holds
            # the individual categorical columns, which have no mean
            agg_result = df.groupby(dimension_name, observed=True)[list(var_measures)].agg(
                [mean, percent_expressed]
            )

            values = []
            dotplot_result = {
                "categories": agg_result.index,
                "name": dimension_name,
                "values": values,
            }
            for var_measure in var_measures:
                series = agg_result[var_measure]
                is_sparse = hasattr(series, "sparse")
                if is_sparse:
                    series = series.sparse.to_dense()
                values.append(
                    {
                        "name": var_measure,
                        "percentExpressed": series["percent_expressed"],
                        "mean": series["mean"],
                    }
                )
            results.append(dotplot_result)
        return results
