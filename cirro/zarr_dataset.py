import pandas as pd
import zarr

from cirro.abstract_dataset import AbstractDataset
from cirro.simple_data import SimpleData


class ZarrDataset(AbstractDataset):

    def __init__(self):
        super().__init__()

    def get_suffixes(self):
        return ['zarr']

    def schema(self, file_system, path):
        store = zarr.open(file_system.get_mapper(path), 'r')
        obs = []
        result = {'version': '1.0.0'}
        obs_keys = store['obs'].keys()
        obs_cat = list(store['obs/__categories'].keys())
        for key in obs_keys:
            if key != '__categories' and key not in obs_cat:
                obs.append(key)
        result['var'] = list(store['var'].index[()])
        result['obs'] = obs
        result['obsCat'] = obs_cat
        result['shape'] = store['X'].attrs['shape']
        embeddings = []
        obsm_summary = store.get('obsm_summary')
        if obsm_summary is not None:
            for key in obsm_summary.keys():
                tokens = key.split('_')
                dim = int(len(tokens) - 3)
                # e.g. X_pca_2_500_max
                embedding = dict(name=key, dimensions=dim)
                embeddings.append(embedding)
        else:
            obsm_keys = store['obsm']
            for key in obsm_keys:
                if key.startswith('X_'):
                    m = store['obsm'][key]
                    shape = m.shape
                    dim = min(3, shape[1])
                    embedding = m.attrs.asdict()
                    embedding['dimensions'] = dim
                    embeddings.append(embedding)
                else:
                    print('Skipping {}'.format(key))
            result['embeddings'] = embeddings
        return result

    def has_precomputed_stats(self, file_system, path, dataset):
        store = zarr.open(file_system.get_mapper(path), 'r')
        return store.get('stats') is not None

    def read_precomputed_stats(self, file_system, path, obs_keys=[], var_keys=[]):
        store = zarr.open(file_system.get_mapper(path), 'r')
        categories = store['obs/__categories']
        result = {}
        if len(obs_keys) > 0:
            for key in obs_keys:
                node = store['stats/obs/' + key]
                if isinstance(node, zarr.hierarchy.Group):  # continuous
                    result[key] = node.attrs.asdict()
                else:
                    counts = node[...]
                    categories_dset = categories[key]
                    categories = categories_dset[...]
                    result[key] = {'categories': categories, 'counts': counts}
        if len(var_keys) > 0:
            for key in var_keys:
                node = store['stats/X/' + key]  # min, max, etc.
                result[key] = node.attrs.asdict()
        return result

    def read_precomputed_grouped_stats(self, file_system, path, obs_keys=[], var_keys=[]):
        store = zarr.open(file_system.get_mapper(path), 'r')
        categories = store['obs/__categories']
        # grouped_stats/obs/node/X_mean, X_num_expressed
        results = []
        if len(obs_keys) > 0 and len(var_keys) > 0:
            var_index = pd.Index(store['var'].index[()])
            indices = var_index.get_indexer_for(var_keys)
            for key in obs_keys:
                values = []
                X_mean = store['grouped_stats/obs/{}/X_mean'.format(key)]
                X_mean = X_mean.get_orthogonal_selection((slice(None), indices))
                X_fraction_expressed = store['grouped_stats/obs/{}/X_fraction_expressed'.format(key)]
                X_fraction_expressed = X_fraction_expressed.get_orthogonal_selection((slice(None), indices))
                categories_dset = categories[key]
                categories = categories_dset[...]
                result = {'categories': categories, 'name': key, 'values': values}
                results.append(result)
                for i in range(len(var_keys)):
                    values.append({'name': var_keys[i],
                                   'fractionExpressed': X_fraction_expressed[:, i].tolist(),
                                   'mean': X_mean[:, i].tolist()})
        return results

    def read_precomputed_basis(self, file_system, path, obs_keys=[], var_keys=[], basis=None):
        store = zarr.open(file_system.get_mapper(path), 'r')
        # categories = store['obs/__categories']
        result = {'values': {}, 'coordinates': {}}
        basis_group = store['obsm_summary/' + basis['full_name']]
        coords_group = basis_group['coords']
        result['bins'] = coords_group['index'][...]
        coords_array = coords_group['value'][...]  # TODO, extract column of interest only shape=(nbins, ndim)
        for i in range(len(basis['coordinate_columns'])):
            result['coordinates'][basis['coordinate_columns'][i]] = coords_array[:, i]
        if len(var_keys) > 0:
            var_index = pd.Index(store['var'].index[()])
            indices = var_index.get_indexer_for(var_keys)
            X = basis_group['X']
            X = X.get_orthogonal_selection((slice(None), indices))
            for i in range(len(indices)):
                result['values'][var_keys[i]] = X[:, i]

        if len(obs_keys) > 0:
            obs_group = basis_group['obs']
            for obs_key in obs_keys:
                obs_node = obs_group[obs_key]
                if isinstance(obs_node, zarr.hierarchy.Group):  # categorical
                    category_mode = obs_node['mode'][...]
                    category_purity = obs_node['purity'][...]
                    # categories_dset = categories[obs_key]
                    # categories = categories_dset[...]
                    # ordered = categories_dset.attrs.get("ordered", False)
                    # value = pd.Categorical.from_codes(category_mode, categories, ordered=ordered)
                    result['values'][obs_key] = dict(value=category_mode.tolist(), purity=category_purity.tolist())
                else:
                    result['values'][obs_key] = obs_node[...]

        return result

    # read zarr dataset stored in anndata format
    def read(self, file_system, path, obs_keys=[], var_keys=[], basis=None, dataset=None):
        store = zarr.open(file_system.get_mapper(path), 'r')
        X = None
        if len(var_keys) > 0:
            var_index = pd.Index(store['var'].index[()])
            indices = var_index.get_indexer_for(var_keys)
            var_index = var_index[indices]  # include selected features only
            group = store['X']
            encoding = group.attrs['encoding-type']
            if encoding != 'csc_matrix':
                raise ValueError('Unsupported encoding')
            shape = group.attrs['shape']
            import scipy.sparse
            data = group["data"]
            mtx = scipy.sparse.csc_matrix(shape, dtype=data.dtype)
            mtx.data = data
            mtx.indices = group["indices"]
            mtx.indptr = group["indptr"]
            X = mtx[:, indices]
        else:
            var_index = pd.Index([])
        obs = pd.DataFrame()
        if len(obs_keys) > 0:
            categories = store['obs/__categories']
            obs_g = store['obs']
            for key in obs_keys:
                dataset = obs_g[key]
                if isinstance(dataset.attrs.get("categories", None), str):
                    categories_dset = categories[key]
                    categories = categories_dset[...]
                    ordered = categories_dset.attrs.get("ordered", False)
                    obs[key] = pd.Categorical.from_codes(dataset[...], categories, ordered=ordered)
                else:
                    obs[key] = dataset[...]
        if basis is not None and len(basis) > 0:
            for b in basis:
                m = store['obsm'][b['name']]
                columns_to_fetch = b['coordinate_columns']
                for i in range(len(columns_to_fetch)):
                    obs[columns_to_fetch[i]] = m[:, i]
        return SimpleData(X, obs, pd.DataFrame(index=var_index))
