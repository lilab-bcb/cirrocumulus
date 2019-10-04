import anndata
import pandas as pd
import scipy


class H5ADBackend:

    def __init__(self, backed='r'):
        self.path_to_data = {}
        self.backed = backed
        # only works with local files

    def get_df(self, filesystem, path, keys, layout_key=None):
        adata = self.path_to_data.get(path, None)
        if adata is None:
            adata = anndata.read(path, backed=self.backed)
            self.path_to_data[path] = adata
        return self.__get_df(adata, keys, layout_key)

    def schema(self, filesystem, path):
        adata = self.path_to_data.get(path, None)
        if adata is None:
            adata = anndata.read(path, backed=self.backed)
            self.path_to_data[path] = adata
        obs_cat = []
        obs = []
        result = {'version': '1'}
        for key in adata.obs_keys():
            if pd.api.types.is_categorical_dtype(adata.obs[key]):
                obs_cat.append(key)
            else:
                obs.append(key)
        result['var'] = adata.var_names
        result['obs'] = obs
        result['obs_cat'] = obs_cat
        layouts = []
        for key in adata.obsm_keys():
            if key.startswith('X_'):
                layouts.append(dict(name=key, dimensions=adata.obsm[key].shape[1]))
        result['layouts'] = layouts
        return result

    def __get_df(self, adata, keys, layout_key=None):
        is_obs = True
        df = pd.DataFrame(index=adata.obs.index.values if is_obs else adata.var.index.values)
        for i in range(len(keys)):
            key = keys[i]
            if key in adata.var_names and is_obs:
                X = adata.obs_vector(key)
                if scipy.sparse.issparse(X):
                    X = X.toarray()
                values = X
            elif key in adata.obs and is_obs:
                values = adata.obs[key].values
            else:
                raise ValueError('{} not found'.format(key))
            df[key] = values
        if layout_key is not None:
            layout_name = layout_key['name']
            layout_data = adata.obsm[layout_name]
            dimensions = layout_key['dimensions']
            for i in range(dimensions):
                df[layout_name + '_' + str(i + 1)] = layout_data[:, i]

        return df
