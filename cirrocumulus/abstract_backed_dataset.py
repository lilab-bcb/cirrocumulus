import os
from abc import abstractmethod

import pandas as pd
import scipy.sparse
from anndata import AnnData

from cirrocumulus.abstract_dataset import AbstractDataset
from cirrocumulus.sparse_dataset import SparseDataset


# string_dtype = h5py.check_string_dtype(dataset.dtype)
# if (string_dtype is not None) and (string_dtype.encoding == "utf-8"):
#     dataset = dataset.asstr()

class AbstractBackedDataset(AbstractDataset):

    def __init__(self):
        super().__init__()

    @abstractmethod
    def is_group(self, node):
        pass

    @abstractmethod
    def open_group(self, filesystem, path):
        pass

    @abstractmethod
    def slice_dense_array(self, X, indices):
        pass

    def get_schema(self, filesystem, path):
        return super().get_schema(filesystem, os.path.join(path, 'index.json.gz'))

    def get_dataset_info(self, filesystem, path):
        d = {}
        with self.open_group(filesystem, path) as group:
            var_ids = group['var/index'][...]
            if pd.api.types.is_object_dtype(var_ids):
                var_ids = var_ids.astype(str)
            d['var'] = pd.Index(var_ids)
            X = group['X']
            d['shape'] = X.attrs['shape'] if self.is_group(X) else X.shape
            if 'uns' in group and 'modules' in group['uns']:
                module_ids = group['uns/modules/var/id'][...]
                if pd.api.types.is_object_dtype(module_ids):
                    module_ids = module_ids.astype(str)
                d['modules'] = pd.Index(module_ids)
        return d

    def read_dataset(self, filesystem, path, keys=None, dataset=None):
        keys = keys.copy()
        X_keys = keys.pop('X', [])
        obs_keys = keys.pop('obs', [])
        basis_keys = keys.pop('basis', [])
        module_keys = keys.pop('module', [])
        X = None
        obs = None
        var = None
        obsm = {}
        adata_modules = None
        dataset_info = self.get_dataset_info(filesystem, path) if len(X_keys) > 0 or len(module_keys) > 0 else None
        root = self.open_group(filesystem, path)
        if len(X_keys) > 0:
            var_ids = dataset_info['var']
            X_node = root['X']
            if len(X_keys) == 1 and isinstance(X_keys[0], slice):  # special case if slice specified for performance
                get_item = X_keys[0]
                X_keys = var_ids[get_item]
            else:
                get_item = var_ids.get_indexer_for(X_keys)
            if self.is_group(X_node):
                X_node = SparseDataset(X_node)  # sparse
                X = X_node[:, get_item]
            else:  # dense
                X = self.slice_dense_array(X_node, get_item)
            var = pd.DataFrame(index=X_keys)
        if len(module_keys) > 0:
            # stored as dense in modules/X, modules/var
            module_ids = dataset_info['modules']
            module_X_node = root['uns/modules/X']
            if len(module_keys) == 1 and isinstance(module_keys[0],
                                                    slice):  # special case if slice specified for performance
                get_item = module_keys[0]
                module_keys = module_ids[get_item]
            else:
                get_item = module_ids.get_indexer_for(module_keys)
            module_X = self.slice_dense_array(module_X_node, get_item)
            adata_modules = AnnData(X=module_X, var=pd.DataFrame(index=module_keys))
        if len(obs_keys) > 0:
            obs = pd.DataFrame()
            group = root['obs']
            for key in obs_keys:
                if key == 'index':
                    values = group['index'][...]
                    if pd.api.types.is_object_dtype(values):
                        values = values.astype(str)
                else:
                    dataset = group[key]
                    values = dataset[...]
                    if "categories" in dataset.attrs:
                        categories = dataset.attrs["categories"]
                        categories_dset = group[categories]
                        categories = categories_dset[...]
                        if pd.api.types.is_object_dtype(categories):
                            categories = categories.astype(str)
                        ordered = categories_dset.attrs.get("ordered", False)
                        values = pd.Categorical.from_codes(values, categories, ordered=ordered)
                obs[key] = values
        if len(basis_keys) > 0:
            group = root['obsm']
            for key in basis_keys:
                embedding_data = group[key][...]
                obsm[key] = embedding_data
                if X is None:
                    X = scipy.sparse.coo_matrix(([], ([], [])), shape=(embedding_data.shape[0], 0))
        adata = AnnData(X=X, obs=obs, var=var, obsm=obsm)
        if adata_modules is not None:
            adata.uns['X_module'] = adata_modules
        return adata
