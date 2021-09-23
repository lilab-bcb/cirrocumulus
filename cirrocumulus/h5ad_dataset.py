import os

import h5py
import numpy as np

from cirrocumulus.abstract_backed_dataset import AbstractBackedDataset


class H5ADDataset(AbstractBackedDataset):

    def __init__(self):
        super().__init__(['h5adx'])

    def is_group(self, node):
        return isinstance(node, h5py.Group)

    def open_group(self, filesystem, path):
        return h5py.File(filesystem.open(os.path.join(path, 'data.h5ad'), 'rb'), mode='r')

    def slice_dense_array(self, X, indices):
        # indexing elements must be in increasing order
        order = np.argsort(indices)
        ordered = indices[order]
        rev_order = np.argsort(order)
        value = X[:, ordered]

        return value[:, rev_order]
