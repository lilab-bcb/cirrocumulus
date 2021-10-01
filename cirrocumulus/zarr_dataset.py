import json

import zarr

from cirrocumulus.abstract_backed_dataset import AbstractBackedDataset


class ZarrDataset(AbstractBackedDataset):

    def __init__(self):
        super().__init__()

    def get_suffixes(self):
        return ['zarr']

    def is_group(self, node):
        return isinstance(node, zarr.hierarchy.Group)

    def open_group(self, filesystem, path):
        return zarr.open_group(filesystem.get_mapper(path), mode='r')

    def slice_dense_array(self, X, indices):
        return X.get_orthogonal_selection((slice(None), indices))

    def get_schema(self, filesystem, path):
        with zarr.open_group(filesystem.get_mapper(path), mode='r') as g:
            return json.loads(str(g['uns']['cirro-schema'][...]))
