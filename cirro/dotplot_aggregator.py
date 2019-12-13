import numpy as np
import pandas as pd
import scipy.sparse
from natsort import natsorted

from cirro.simple_data import SimpleData


class DotPlotAggregator:

    def __init__(self, var_measures, dimensions):
        self.var_measures = var_measures
        self.dimensions = dimensions

    def execute(self, adata):
        results = []
        # {categories:[], name:'', values:[{name:'', fractionExpressed:0, mean:0}]}

        var_measures = self.var_measures
        dimensions = self.dimensions
        X = adata.X[:, SimpleData.get_var_indices(adata, var_measures)]
        issparse = scipy.sparse.issparse(X)

        for dimension in dimensions:
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

            index = pd.Index(group_names)
            sorted_categories = natsorted(index)
            reordered_indices = index.get_indexer_for(sorted_categories)
            values = []
            dotplot_result = {'categories': sorted_categories, 'name': dimension, 'values': values}
            mean_output = mean_output[reordered_indices]
            fraction_expressed_output = fraction_expressed_output[reordered_indices]
            for i in range(mean_output.shape[1]):
                name = var_measures[i]
                values.append({'name': name,
                               'fractionExpressed': fraction_expressed_output[:, i].tolist(),
                               'mean': mean_output[:, i].tolist()})
            results.append(dotplot_result)
        return results
