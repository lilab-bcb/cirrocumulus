import json

import zarr

from cirrocumulus.abstract_backed_dataset import AbstractBackedDataset
from cirrocumulus.anndata_util import dataset_schema


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
        g = zarr.open_group(filesystem.get_mapper(path), mode='r')
        if 'cirro-schema' in g['uns']:
            s = str(g['uns']['cirro-schema'][()])
            # hack to fix encoding error in old version
            s = s.replace('"id"›', '"id":')
            return json.loads(s)
        return dataset_schema(g, n_features=0)
