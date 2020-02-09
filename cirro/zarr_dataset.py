import pandas as pd
import zarr

from cirro.abstract_dataset import AbstractDataset


class ZarrDataset(AbstractDataset):

    def __init__(self):
        super().__init__()

    def get_suffixes(self):
        return ['zarr', 'zjson']

    def to_pandas(self, file_system, path, columns=None):
        store = file_system.get_mapper(path)
        g = zarr.open(store)
        result = {}
        keys = columns if columns is not None else g.array_keys(False)
        for key in keys:
            result[key] = g[key][()]
        return pd.DataFrame.from_dict(result)
