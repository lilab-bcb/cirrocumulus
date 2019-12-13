import numpy as np
import scipy.sparse

from cirro.simple_data import SimpleData


def get_basis(basis):
    if basis is not None:
        embedding_ndim = 2
        if basis.endswith('3d'):
            basis = basis[0:len(basis) - 3]
            embedding_ndim = 3
        coordinate_columns = []
        for i in range(embedding_ndim):
            coordinate_columns.append(basis + '_' + str(i + 1))
        return {'name': basis, 'dimensions': embedding_ndim, 'coordinate_columns': coordinate_columns}


class EmbeddingAggregator:

    def __init__(self, obs_measures, var_measures, dimensions, count, nbins, basis, agg_function):
        self.nbins = nbins
        self.agg_function = agg_function
        self.obs_measures = obs_measures
        self.var_measures = var_measures
        self.dimensions = dimensions
        self.count = count
        self.add_count = count or nbins is not None
        self.basis = basis

    @staticmethod
    def get_bin_level_agg_dict(measures, coordinate_columns, agg_function):
        agg_func = {}
        for column in measures:
            if column == '__count':
                agg_func[column] = 'sum'
            else:
                agg_func[column] = agg_function
        for column in coordinate_columns:
            agg_func[column] = 'min'
        return agg_func


    @staticmethod
    def convert_coords_to_bin(df, nbins, coordinate_columns, coordinate_column_to_range=None):
        # replace coordinates with bin, set df index to bin
        for name in coordinate_columns:
            values = df[name].values
            if coordinate_column_to_range is not None:
                view_column_range = coordinate_column_to_range[name]
                column_min = view_column_range[0]
                column_max = view_column_range[1]
            else:
                column_min = values.min()
                column_max = values.max()
            df[name] = np.floor(np.interp(values, [column_min, column_max], [0, nbins - 1])).astype(int)
        if len(coordinate_columns) == 2:
            df.index = df[coordinate_columns[0]] * nbins + df[coordinate_columns[1]]
        else:
            df.index = df[coordinate_columns[2]] + nbins * (
                    df[coordinate_columns[1]] + nbins * df[coordinate_columns[0]])

    def execute(self, adata):
        result = {'coordinates': {}, 'values': {}}
        var_measures = self.var_measures
        obs_measures = self.obs_measures
        dimensions = self.dimensions
        add_count = self.add_count
        basis = self.basis

        if self.nbins is not None:
            agg_function = self.agg_function
            df = adata.obs
            if add_count:
                df['__count'] = 1.0
            # bin level summary, coordinates have already been converted
            obs_measure_agg_dict = {}
            for column in basis['coordinate_columns']:
                obs_measure_agg_dict[column] = 'min'
            for column in obs_measures:
                obs_measure_agg_dict[column] = agg_function
            grouped = df.groupby(df.index)
            X_output = None
            has_var_measures = len(var_measures) > 0
            has_dimensions = len(dimensions) > 0
            if has_var_measures:
                X = adata.X[:, SimpleData.get_var_indices(adata, var_measures)]
                is_sparse = scipy.sparse.issparse(X)
            obs_summary = grouped.agg(obs_measure_agg_dict)
            dimension_purity_output = {}
            dimension_mode_output = {}
            for column in dimensions:
                dimension_purity_output[column] = []
                dimension_mode_output[column] = []
            if has_var_measures or len(dimensions) > 0:
                for key, g in grouped:
                    indices = grouped.indices[key]
                    if has_dimensions:
                        group_df = df.iloc[indices]
                        for dimension in dimensions:
                            value_counts = group_df[dimension].value_counts(sort=False)
                            largest = value_counts.nlargest(1)
                            purity = largest[0] / value_counts.sum()
                            dimension_mode_output[dimension].append(largest.index[0])
                            dimension_purity_output[dimension].append(purity)
                    if has_var_measures:
                        X_group = X[indices]
                        if agg_function == 'max':
                            X_summary = X_group.max(axis=0)
                            if is_sparse:
                                X_summary = X_summary.toarray().flatten()
                        elif agg_function == 'min':
                            X_summary = X_group.min(axis=0)
                            if is_sparse:
                                X_summary = X_summary.toarray().flatten()
                        elif agg_function == 'mean':
                            X_summary = X_group.mean(axis=0)
                            if is_sparse:
                                X_summary = X_summary.A1
                        elif agg_function == 'sum':
                            X_summary = X_group.sum(axis=0)
                            if is_sparse:
                                X_summary = X_summary.A1
                        X_output = np.vstack((X_output, X_summary)) if X_output is not None else X_summary
                for i in range(len(var_measures)):
                    result['values'][var_measures[i]] = X_output[:, i].tolist()
            for i in range(len(obs_measures)):
                result['values'][obs_measures[i]] = obs_summary[obs_measures[i]].tolist()
            for column in dimensions:
                result['values'][column] = dict(mode=dimension_mode_output[column],
                    purity=dimension_purity_output[column])
            result['bins'] = obs_summary.index.to_list()
            for column in basis['coordinate_columns']:
                result['coordinates'][column] = obs_summary[column].tolist()
        else:
            if add_count:
                result['values']['__count'] = np.ones(adata.shape[0]).tolist()
            if len(var_measures) > 0:
                X = adata.X[:, SimpleData.get_var_indices(adata, var_measures)]
                is_sparse = scipy.sparse.issparse(X)
                for i in range(len(var_measures)):
                    result['values'][var_measures[i]] = X[:, i].tolist() if not is_sparse else X[:,
                                                                                               i].toarray().flatten().tolist()
            for column in obs_measures + dimensions:
                result['values'][column] = adata.obs[column].to_list()
            for column in basis['coordinate_columns']:
                result['coordinates'][column] = adata.obs[column].to_list()
        return result
