from collections.abc import Mapping
from enum import Enum
from functools import singledispatch, _find_impl
from types import MappingProxyType
from typing import Callable, Type, TypeVar
from warnings import warn

import numcodecs
import numpy as np
import pandas as pd
import zarr
from packaging import version
from pandas.api.types import is_categorical_dtype
from scipy import sparse


class WriteWarning(UserWarning):
    pass


T = TypeVar("T")


def _to_fixed_length_strings(value: np.ndarray) -> np.ndarray:
    """\
    Convert variable length strings to fixed length.

    Currently a workaround for
    https://github.com/zarr-developers/zarr-python/pull/422
    """
    new_dtype = []
    for dt_name, (dt_type, dt_offset) in value.dtype.fields.items():
        if dt_type.kind == "O":
            #  Assuming the objects are str
            size = max(len(x.encode()) for x in value.getfield("O", dt_offset))
            new_dtype.append((dt_name, ("U", size)))
        else:
            new_dtype.append((dt_name, dt_type))
    return value.astype(new_dtype)


def check_key(key):
    """Checks that passed value is a valid h5py key.

    Should convert it if there is an obvious conversion path, error otherwise.
    """
    typ = type(key)
    if issubclass(typ, str):
        return str(key)
    # TODO: Should I try to decode bytes? It's what h5py would do,
    # but it will be read out as a str.
    # elif issubclass(typ, bytes):
    # return key
    else:
        raise TypeError(f"{key} of type {typ} is an invalid key. Should be str.")


# -------------------------------------------------------------------------------
# Generic functions
# -------------------------------------------------------------------------------


@singledispatch
def write_attribute(*args, **kwargs):
    raise NotImplementedError("Unrecognized argument types for `write_attribute`.")


class EncodingVersions(Enum):
    raw = "0.1.0"
    csr_matrix = csc_matrix = "0.1.0"
    dataframe = "0.1.0"

    def check(self, key: str, encoded_version: str):
        if version.parse(encoded_version) > version.parse(self.value):
            warn(
                f"The supported version for decoding {self.name} is {self.value}, "
                f"but a {self.name} with version {encoded_version} "
                f"was encountered at {key}.",
                FutureWarning,
            )


def _write_method(cls: Type[T]) -> Callable[[zarr.Group, str, T], None]:
    return _find_impl(cls, ZARR_WRITE_REGISTRY)


@write_attribute.register(zarr.Group)
def write_attribute_zarr(f, key, value, dataset_kwargs=MappingProxyType({})):
    if key in f:
        del f[key]
    _write_method(type(value))(f, key, value, dataset_kwargs)


def write_mapping(f, key, value: Mapping, dataset_kwargs=MappingProxyType({})):
    for sub_k, sub_v in value.items():
        if not isinstance(key, str):
            warn(
                f"dict key {key} transformed to str upon writing to zarr, using "
                "string keys is recommended.",
                WriteWarning,
            )
        write_attribute(f, f"{key}/{sub_k}", sub_v, dataset_kwargs)


def write_dataframe(z, key, df, dataset_kwargs=MappingProxyType({})):
    # Check arguments
    for reserved in ("__categories", "_index"):
        if reserved in df.columns:
            raise ValueError(f"{reserved!r} is a reserved name for dataframe columns.")

    col_names = [check_key(c) for c in df.columns]

    if df.index.name is not None:
        index_name = df.index.name
    else:
        index_name = "_index"
    index_name = check_key(index_name)

    group = z.create_group(key)
    group.attrs["encoding-type"] = "dataframe"
    group.attrs["encoding-version"] = EncodingVersions.dataframe.value
    group.attrs["column-order"] = col_names
    group.attrs["_index"] = index_name

    write_series(group, index_name, df.index, dataset_kwargs)
    for col_name, (_, series) in zip(col_names, df.items()):
        write_series(group, col_name, series, dataset_kwargs)


def write_series(group, key, series, dataset_kwargs=MappingProxyType({})):
    if series.dtype == object:
        group.create_dataset(
            key,
            shape=series.shape,
            dtype=object,
            object_codec=numcodecs.VLenUTF8(),
            **dataset_kwargs,
        )
        group[key][:] = series.values
    elif is_categorical_dtype(series):
        # This should work for categorical Index and Series
        categorical: pd.Categorical = series.values
        categories: np.ndarray = categorical.categories.values
        codes: np.ndarray = categorical.codes
        category_key = f"__categories/{key}"

        write_array(group, category_key, categories, dataset_kwargs=dataset_kwargs)
        write_array(group, key, codes, dataset_kwargs=dataset_kwargs)

        group[key].attrs["categories"] = category_key
        # Must coerce np.bool_ to bool for json writing
        group[category_key].attrs["ordered"] = bool(categorical.ordered)
    else:
        write_array(group, key, series.values, dataset_kwargs=dataset_kwargs)


def write_not_implemented(f, key, value, dataset_kwargs=MappingProxyType({})):
    # If it’s not an array, try and make it an array. If that fails, pickle it.
    # Maybe rethink that, maybe this should just pickle,
    # and have explicit implementations for everything else
    raise NotImplementedError(
        f"Failed to write value for {key}, since a writer for type {type(value)}"
        f" has not been implemented yet."
    )


def write_list(g, key, value, dataset_kwargs=MappingProxyType({})):
    write_array(g, key, np.array(value), dataset_kwargs)


def write_array(g, key, value, dataset_kwargs=MappingProxyType({})):
    if value.dtype == object:
        g.create_dataset(
            key,
            shape=value.shape,
            dtype=object,
            object_codec=numcodecs.VLenUTF8(),
            **dataset_kwargs,
        )
        g[key][:] = value
    elif value.dtype.kind == "V":
        # Structured dtype
        g.create_dataset(key, data=_to_fixed_length_strings(value), **dataset_kwargs)
    else:
        g.create_dataset(key, data=value, **dataset_kwargs)


# TODO: Not working quite right

def write_scalar(f, key, value, dataset_kwargs=MappingProxyType({})):
    f.create_dataset(key, data=np.array(value), **dataset_kwargs)


def write_none(f, key, value, dataset_kwargs=MappingProxyType({})):
    pass


# TODO: Figure out what to do with dataset_kwargs for these

def write_csr(f, key, value: sparse.csr_matrix, dataset_kwargs=MappingProxyType({})):
    group = f.create_group(key)
    group.attrs["encoding-type"] = "csr_matrix"
    group.attrs["encoding-version"] = EncodingVersions.csr_matrix.value
    group.attrs["shape"] = value.shape
    write_array(group, "data", value.data, dataset_kwargs=dataset_kwargs)
    write_array(group, "indices", value.indices, dataset_kwargs=dataset_kwargs)
    write_array(group, "indptr", value.indptr, dataset_kwargs=dataset_kwargs)


def write_csc(f, key, value: sparse.csc_matrix, dataset_kwargs=MappingProxyType({})):
    group = f.create_group(key)
    group.attrs["encoding-type"] = "csc_matrix"
    group.attrs["encoding-version"] = EncodingVersions.csc_matrix.value
    group.attrs["shape"] = value.shape
    write_array(group, "data", value.data, dataset_kwargs=dataset_kwargs)
    write_array(group, "indices", value.indices, dataset_kwargs=dataset_kwargs)
    write_array(group, "indptr", value.indptr, dataset_kwargs=dataset_kwargs)


def write_raw(f, key, value, dataset_kwargs=MappingProxyType({})):
    group = f.create_group(key)
    group.attrs["encoding-type"] = "raw"
    group.attrs["encoding-version"] = EncodingVersions.raw.value
    group.attrs["shape"] = value.shape
    write_attribute(group, "X", value.X, dataset_kwargs)
    write_attribute(group, "var", value.var, dataset_kwargs)
    write_attribute(group, "varm", value.varm, dataset_kwargs)


ZARR_WRITE_REGISTRY = {
    type(None): write_none,
    Mapping: write_mapping,
    object: write_not_implemented,
    np.ndarray: write_array,  # Possibly merge with write_series
    list: write_list,
    pd.DataFrame: write_dataframe,
    # Raw: write_raw,
    # object: write_not_implemented,
    # h5py.Dataset: write_basic,
    # type(None): write_none,
    str: write_scalar,
    float: write_scalar,
    np.floating: write_scalar,
    bool: write_scalar,
    np.bool_: write_scalar,
    int: write_scalar,
    np.integer: write_scalar,
    sparse.csr_matrix: write_csr,
    sparse.csc_matrix: write_csc,
}
