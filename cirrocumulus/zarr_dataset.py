import json

import zarr

from cirrocumulus.abstract_backed_dataset import AbstractBackedDataset
from cirrocumulus.anndata_util import dataset_schema


class AttributeGroup:
    """Expose a zarr group's members as attributes.

    zarr 2 groups did this themselves and zarr 3 dropped it, but ``dataset_schema`` is
    written in that style because it also accepts a plain ``AnnData``.
    """

    def __init__(self, group):
        self._group = group

    def __getattr__(self, name):
        try:
            return self._group[name]
        except KeyError as err:
            raise AttributeError(name) from err

    def __getitem__(self, key):
        return self._group[key]

    def __contains__(self, key):
        return key in self._group


class ZarrDataset(AbstractBackedDataset):
    def __init__(self):
        super().__init__()

    def get_suffixes(self):
        return ["zarr"]

    def is_group(self, node):
        return isinstance(node, zarr.Group)

    def open_group(self, filesystem, path):
        return zarr.open_group(filesystem.get_mapper(path), mode="r")

    def slice_dense_array(self, X, indices):
        return X.get_orthogonal_selection((slice(None), indices))

    def get_schema(self, filesystem, path):
        g = self.open_group(filesystem, path)
        if "uns" in g and "cirro-schema" in g["uns"]:
            return json.loads(str(g["uns"]["cirro-schema"][()]))
        return dataset_schema(AttributeGroup(g), n_features=0)
