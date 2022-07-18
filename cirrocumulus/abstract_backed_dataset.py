from abc import abstractmethod

import pandas as pd
import scipy.sparse
from anndata import AnnData

from cirrocumulus.abstract_dataset import AbstractDataset
from cirrocumulus.anndata_util import ADATA_LAYERS_UNS_KEY, ADATA_MODULE_UNS_KEY
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

    def get_result(self, filesystem, path, dataset, result_id):
        g = self.open_group(filesystem, path)
        uns = g["uns"]
        if result_id in uns:
            return str(uns[result_id][...])
        return super().get_result(filesystem, path, dataset, result_id)

    def get_dataset_info(self, filesystem, path):
        d = {}
        root = self.open_group(filesystem, path)
        var_group = root["var"]
        var_group_index_field = var_group.attrs["_index"]
        var_ids = var_group[var_group_index_field][...]
        if pd.api.types.is_object_dtype(var_ids):
            var_ids = var_ids.astype(str)
        d["var"] = pd.Index(var_ids)
        X = root["X"]
        d["shape"] = X.attrs["shape"] if self.is_group(X) else X.shape
        if "layers" in root:
            d["layers"] = root["layers"].keys()
        if "uns" in root:
            uns_group = root["uns"]
            if "module" in uns_group:
                module_var_group = uns_group["module/var"]
                module_var_group_index_field = module_var_group.attrs["_index"]
                module_ids = module_var_group[module_var_group_index_field][...]
                if pd.api.types.is_object_dtype(module_ids):
                    module_ids = module_ids.astype(str)
                d["module"] = pd.Index(module_ids)
        return d

    def get_X(self, var_ids, keys, node):
        if len(keys) == 1 and isinstance(
            keys[0], slice
        ):  # special case if slice specified for performance
            get_item = keys[0]
            keys = var_ids[get_item]
        else:
            get_item = var_ids.get_indexer_for(keys)

        if self.is_group(node):
            sparse_dataset = SparseDataset(node)  # sparse
            X = sparse_dataset[:, get_item]
        else:  # dense
            X = self.slice_dense_array(node, get_item)
        var = pd.DataFrame(index=keys)
        return X, var

    def read_dataset(self, filesystem, path, keys=None, dataset=None):
        keys = keys.copy()
        X_keys = keys.pop("X", [])
        obs_keys = keys.pop("obs", [])
        basis_keys = keys.pop("basis", [])
        module_keys = keys.pop("module", [])
        # additional keys belong to layers
        X = None
        obs = None
        var = None
        obsm = {}
        adata_modules = None
        dataset_info = self.get_dataset_info(filesystem, path)
        root = self.open_group(filesystem, path)
        layers = {}
        for layer_key in keys.keys():
            X_layer, var_layer = self.get_X(
                dataset_info["var"], keys[layer_key], root["layers"][layer_key]
            )
            adata_layer = AnnData(X=X_layer, var=var_layer)
            layers[layer_key] = adata_layer
        if len(X_keys) > 0:
            X, var = self.get_X(dataset_info["var"], X_keys, root["X"])
        if len(obs_keys) > 0:
            obs = pd.DataFrame(index=pd.RangeIndex(dataset_info["shape"][0]).astype(str))
            group = root["obs"]
            for key in obs_keys:
                if key == "index":
                    index_field = group.attrs["_index"]
                    values = group[index_field][...]
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
        if len(module_keys) > 0:
            module_ids = dataset_info["module"]
            module_X_node = root["uns/module/X"]
            module_X, module_var = self.get_X(module_ids, module_keys, module_X_node)
            adata_modules = AnnData(X=module_X, var=module_var, obs=obs)  # obs is shared
        if len(basis_keys) > 0:
            group = root["obsm"]
            for key in basis_keys:
                embedding_data = group[key][...]
                obsm[key] = embedding_data
                if X is None:
                    X = scipy.sparse.coo_matrix(([], ([], [])), shape=(embedding_data.shape[0], 0))
        if X is None and obs is None and len(obsm.keys()) == 0:
            if dataset_info is None:
                dataset_info = self.get_dataset_info(filesystem, path)
            obs = pd.DataFrame(index=pd.RangeIndex(dataset_info["shape"][0]).astype(str))
        adata = AnnData(X=X, obs=obs, var=var, obsm=obsm)
        if adata_modules is not None:
            adata.uns[ADATA_MODULE_UNS_KEY] = adata_modules
        adata.uns[ADATA_LAYERS_UNS_KEY] = layers
        return adata
