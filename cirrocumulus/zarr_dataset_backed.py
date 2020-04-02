from itertools import accumulate, chain
from typing import Tuple, Sequence, Iterable

import numpy as np
import pandas as pd
import scipy.sparse as ss
import zarr
from scipy.sparse import _sparsetools
from scipy.sparse.compressed import _cs_matrix

from cirrocumulus.abstract_dataset import AbstractDataset
from cirrocumulus.simple_data import SimpleData


class BackedSparseMatrix(_cs_matrix):


    def _set_many(self, i: Iterable[int], j: Iterable[int], x):
        """\
        Sets value at each (i, j) to x
        Here (i,j) index major and minor respectively,
        and must not contain duplicate entries.
        """
        # Scipy 1.3+ compat
        n_samples = 1 if np.isscalar(x) else len(x)
        offsets = self._offsets(i, j, n_samples)

        if -1 not in offsets:
            # make a list for interaction with h5py
            offsets = list(offsets)
            # only affects existing non-zero cells
            self.data[offsets] = x
            return

        else:
            raise ValueError(
                "You cannot change the sparsity structure of a SparseDataset."
            )
            # replace where possible
            # mask = offsets > -1
            # # offsets[mask]
            # bool_data_mask = np.zeros(len(self.data), dtype=bool)
            # bool_data_mask[offsets[mask]] = True
            # self.data[bool_data_mask] = x[mask]
            # # self.data[offsets[mask]] = x[mask]
            # # only insertions remain
            # mask = ~mask
            # i = i[mask]
            # i[i < 0] += M
            # j = j[mask]
            # j[j < 0] += N
            # self._insert_many(i, j, x[mask])

    def _zero_many(self, i: Sequence[int], j: Sequence[int]):
        """\
        Sets value at each (i, j) to zero, preserving sparsity structure.
        Here (i,j) index major and minor respectively.
        """
        offsets = self._offsets(i, j, len(i))

        # only assign zeros to the existing sparsity structure
        self.data[list(offsets[offsets > -1])] = 0

    def _offsets(
            self, i: Iterable[int], j: Iterable[int], n_samples: int
    ) -> np.ndarray:
        i, j, M, N = self._prepare_indices(i, j)
        offsets = np.empty(n_samples, dtype=self.indices.dtype)
        ret = _sparsetools.csr_sample_offsets(
            M, N, self.indptr, self.indices, n_samples, i, j, offsets
        )
        if ret == 1:
            # rinse and repeat
            self.sum_duplicates()
            _sparsetools.csr_sample_offsets(
                M, N, self.indptr, self.indices, n_samples, i, j, offsets
            )
        return offsets


class backed_csc_matrix(BackedSparseMatrix, ss.csc_matrix):
    def _get_sliceXint(self, row: slice, col: int) -> ss.csc_matrix:
        return ss.csc_matrix(
            get_compressed_vector(self, col), shape=(self.shape[0], 1)
        )[row, :]

    def _get_sliceXslice(self, row: slice, col: slice) -> ss.csc_matrix:
        out_shape = (
                slice_len(row, self.shape[0]),
                slice_len(col, self.shape[1]),
        )
        if out_shape[1] == 1:
            return self._get_sliceXint(row, slice_as_int(col, self.shape[1]))
        elif out_shape[0] == self.shape[0] and out_shape[1] < self.shape[1]:
            return self._get_sliceXarray(row, np.arange(*col.indices(self.shape[1])))
        return super()._get_sliceXslice(row, col)

    def _get_sliceXarray(self, row: slice, col: Sequence[int]) -> ss.csc_matrix:
        idxs = np.asarray(col)
        if idxs.dtype == bool:
            idxs = np.where(idxs)
        return ss.csc_matrix(
            get_compressed_vectors(self, idxs), shape=(self.shape[0], len(idxs))
        )[row, :]


def slice_len(s: slice, l: int) -> int:
    """Returns length of `a[s]` where `len(a) == l`."""
    return len(range(*s.indices(l)))


def slice_as_int(s: slice, l: int) -> int:
    """Converts slices of length 1 to the integer index they’ll access."""
    out = list(range(*s.indices(l)))
    assert len(out) == 1
    return out[0]


def get_compressed_vectors(
        x: BackedSparseMatrix, row_idxs: Iterable[int]
) -> Tuple[Sequence, Sequence, Sequence]:
    slices = [slice(*(x.indptr[i: i + 2])) for i in row_idxs]
    data = np.concatenate([x.data[s] for s in slices])
    indices = np.concatenate([x.indices[s] for s in slices])
    indptr = list(accumulate(chain((0,), (s.stop - s.start for s in slices))))
    return data, indices, indptr


def get_compressed_vector(
        x: BackedSparseMatrix, idx: int
) -> Tuple[Sequence, Sequence, Sequence]:
    s = slice(*(x.indptr[idx: idx + 2]))
    data = x.data[s]
    indices = x.indices[s]
    indptr = [0, len(data)]
    return data, indices, indptr


class ZarrDatasetBacked(AbstractDataset):

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
        result['var'] = store['var/_index'][...].tolist()
        result['obs'] = obs
        result['obsCat'] = obs_cat
        result['shape'] = store['X'].attrs['shape']
        embeddings = []
        obsm_summary = store.get('obsm_summary')
        if obsm_summary is not None:
            for key in obsm_summary.keys():
                embedding = obsm_summary[key].attrs.asdict()
                embeddings.append(embedding)
        else:
            obsm_keys = store['obsm']
            for key in obsm_keys:
                if key.startswith('X_'):
                    m = store['obsm'][key]
                    shape = m.shape
                    dim = min(3, shape[1])
                    embedding = m.attrs.asdict()
                    embedding['name'] = key
                    embedding['dimensions'] = dim
                    embeddings.append(embedding)
                    if shape[1] >= 3:
                        embedding = embedding.copy()
                        embedding['dimensions'] = 2
                        embeddings.append(embedding)
                else:
                    print('Skipping {}'.format(key))
        result['embeddings'] = embeddings
        return result

    def read_precomputed_basis(self, file_system, path, obs_keys=[], var_keys=[], basis=None):
        store = zarr.open(file_system.get_mapper(path), 'r')
        # categories = store['obs/__categories']
        result = {'values': {}, 'coordinates': {}}
        basis_group = store['obsm_summary/' + basis['full_name']]
        coords_group = basis_group['coords']
        result['bins'] = coords_group['index'][...].tolist()
        coords_array = coords_group['value'][...]  # TODO, extract column of interest only shape=(nbins, ndim)
        for i in range(len(basis['coordinate_columns'])):
            result['coordinates'][basis['coordinate_columns'][i]] = coords_array[:, i].tolist()
        if len(var_keys) > 0:
            var_index = pd.Index(store['var/_index'][...])
            indices = var_index.get_indexer_for(var_keys)
            X = basis_group['X']
            X = X.get_orthogonal_selection((slice(None), indices))
            for i in range(len(indices)):
                result['values'][var_keys[i]] = X[:, i].tolist()

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
                    result['values'][obs_key] = obs_node[...].tolist()

        return result

    # read zarr dataset stored in anndata format
    def read(self, file_system, path, obs_keys=[], var_keys=[], basis=None, dataset=None, schema=None):
        store = zarr.open(file_system.get_mapper(path), 'r')
        X = None
        if len(var_keys) > 0:
            var_index = pd.Index(store['var/_index'][...])
            indices = var_index.get_indexer_for(var_keys)
            var_index = var_index[indices]  # include selected features only
            group = store['X']
            encoding = group.attrs['encoding-type']
            if encoding != 'csc_matrix':
                raise ValueError('Unsupported encoding')

            mtx = backed_csc_matrix(tuple(group.attrs['shape']), dtype=group["data"].dtype)
            mtx.data = group["data"]
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
                if b['precomputed']:
                    basis_group = store['obsm_summary/' + b['full_name']]
                    coords_group = basis_group['coords']
                    obs[b['full_name']] = coords_group['bins'][...].tolist()
                    coords_array = coords_group['bin_coords'][
                        ...]  # TODO, extract column of interest only shape=(nbins, ndim)
                    for i in range(len(b['coordinate_columns'])):
                        obs[b['coordinate_columns'][i]] = coords_array[:, i].tolist()
                else:
                    m = store['obsm'][b['name']]
                    columns_to_fetch = b['coordinate_columns']
                    for i in range(len(columns_to_fetch)):
                        obs[columns_to_fetch[i]] = m[:, i]
        return SimpleData(X, obs, pd.DataFrame(index=var_index))
