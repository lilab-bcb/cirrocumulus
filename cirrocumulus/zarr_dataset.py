import zarr

from cirrocumulus.abstract_backed_dataset import AbstractBackedDataset


class ZarrDataset(AbstractBackedDataset):

    def __init__(self):
        super().__init__(['zarr'])

    def is_group(self, node):
        return isinstance(node, zarr.hierarchy.Group)

    def open_group(self, filesystem, path):
        return zarr.open_group(filesystem.get_mapper(path), mode='r')

    def slice_dense_array(self, X, indices):
        return X.get_orthogonal_selection((slice(None), indices))
