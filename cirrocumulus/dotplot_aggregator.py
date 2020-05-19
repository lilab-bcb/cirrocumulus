import numpy as np
import pandas as pd
import scipy.sparse
from cirrocumulus.simple_data import SimpleData
from natsort import natsorted
from pandas import CategoricalDtype


class DotPlotAggregator:

    def __init__(self, var_measures, dimensions):
        self.var_measures = var_measures
        self.dimensions = dimensions

    def execute(self, adata):
        results = []
        # {categories:[], name:'', values:[{name:'', fractionExpressed:0, mean:0}]}
        var_measures = self.var_measures
        dimensions = self.dimensions
        if len(var_measures) == 0 or len(dimensions) == 0:
            return results
        X = adata.X[:, SimpleData.get_var_indices(adata, var_measures)]
        issparse = scipy.sparse.issparse(X)
        df = adata.obs
        for dimension in dimensions:
            if pd.api.types.is_categorical_dtype(df[dimension]):
                if not df[dimension].dtype.ordered:
                    df[dimension] = df[dimension].astype(
                        CategoricalDtype(natsorted(df[dimension].dtype.categories), ordered=True))
                if len(df[dimension].dtype.categories) <= 1:
                    continue
            grouped = adata.obs.groupby(dimension)
            group_names = []
            mean_output = None
            fraction_expressed_output = None
            for key, g in grouped:
                group_names.append(key)
                indices = grouped.indices[key]
                X_group = X[indices]
                mean_values = X_group.mean(axis=0)

                if issparse:
                    mean_values = mean_values.A1
                    fraction_expressed = X_group.getnnz(axis=0) / X_group.shape[0]
                else:
                    fraction_expressed = (X_group != 0).sum(axis=0) / (X_group.shape[0])

                mean_output = np.vstack((mean_output, mean_values)) if mean_output is not None else mean_values
                fraction_expressed_output = np.vstack((fraction_expressed_output,
                                                       fraction_expressed)) if fraction_expressed_output is not None else fraction_expressed

            values = []
            dotplot_result = {'categories': group_names, 'name': dimension, 'values': values}
            for i in range(mean_output.shape[1]):
                name = var_measures[i]
                values.append({'name': name,
                               'fractionExpressed': fraction_expressed_output[:, i],
                               'mean': mean_output[:, i]})
            results.append(dotplot_result)
        return results
