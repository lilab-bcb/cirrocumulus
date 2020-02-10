import logging

import numcodecs
import numpy as np
import pandas as pd
import zarr

from cirro.simple_data import SimpleData

logger = logging.getLogger("cirro")
from cirro.abstract_dataset import AbstractDataset


def write_zarr(d, output, name):
    g = output.create_group(name)
    for key in d:
        value = d[key]
        if isinstance(value, list):
            value = np.array(value)
        dtype = value.dtype
        object_codec = None
        if pd.api.types.is_categorical_dtype(value):
            object_codec = numcodecs.Categorize(value.cat.categories, dtype=object)
            dtype = object
            value = value.astype(str)
        elif dtype == object:
            object_codec = numcodecs.MsgPack()
        if isinstance(value, pd.Series):
            value = value.values
        ds = g.create_dataset(key, shape=value.shape, chunks=None, dtype=dtype, object_codec=object_codec)
        ds[:] = value


class ZarrDataset(AbstractDataset):

    def __init__(self):
        super().__init__()

    # def get_dataset_attrs(self, file_system, path):
    #     with file_system.open(os.path.join(path, 'data', 'index.json')) as s:
    #         return json.load(s)

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
        obsm_keys = store['obsm']
        for key in obsm_keys:
            if key.startswith('X_'):
                shape = store['obsm'][key].shape
                dim = min(3, shape[1])
                embedding = dict(name=key, dimensions=dim)
                if dim == 3:
                    embeddings.append(embedding)
                    embedding = embedding.copy()
                    embedding['dimensions'] = 2
                    embeddings.append(embedding)
                else:
                    embeddings.append(embedding)
            else:
                print('Skipping {}'.format(key))
        result['embeddings'] = embeddings
        return result


    # read zarr dataset stored in anndata format
    def read(self, file_system, path, obs_keys=[], var_keys=[], basis=None, dataset=None):
        store = zarr.open(file_system.get_mapper(path), 'r')
        X = None
        if len(var_keys) > 0:
            var_index = pd.Index(store['var'].index[()])
            indices = var_index.get_indexer_for(var_index)
            var_index = var_index[indices]
            X = store['X'][:, indices][()]
        else:
            var_index = pd.Index([])
        obs = pd.DataFrame()
        if len(obs_keys) > 0:
            categories = store['obs/__categories']
            obs_g = store['obs']
            for key in obs_keys:
                dataset = obs_g[key]
                if isinstance(dataset.attrs["categories"], str):
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


    def to_pandas(self, file_system, path, columns=None):
        store = file_system.get_mapper(path)
        g = zarr.open(store)
        result = {}
        if isinstance(g, zarr.hierarchy.Group):
            keys = columns if columns is not None else g.array_keys(False)
        else:
            keys = columns if columns is not None else g.dtype.fields.keys()
        for key in keys:
            result[key] = g[key][()]
        return pd.DataFrame.from_dict(result)
