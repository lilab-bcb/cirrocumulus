import os

import fsspec
import numpy as np
import pandas as pd
import pytest
import scipy.sparse

from cirrocumulus.parquet_dataset import ParquetDataset
from cirrocumulus.prepare_data import PrepareData
from cirrocumulus.zarr_dataset import ZarrDataset


def read_and_diff(ds_reader, path, test_data, measures, dimensions, continuous_obs, basis):
    dataset = dict(id='')
    fs = fsspec.filesystem('file')
    prepared_adata = ds_reader.read_dataset(filesystem=fs, path=path, dataset=dataset,
                                            keys=dict(X=measures, obs=dimensions + continuous_obs,
                                                      basis=[basis]))

    assert scipy.sparse.issparse(test_data.X) == scipy.sparse.issparse(prepared_adata.X)
    if scipy.sparse.issparse(test_data.X):
        test_data.X = test_data.X.toarray()
        prepared_adata.X = prepared_adata.X.toarray()
    np.testing.assert_equal(test_data.X, prepared_adata.X)
    for key in dimensions + continuous_obs:
        pd.testing.assert_series_equal(test_data.obs[key], prepared_adata.obs[key], check_index=False,
                                       check_flags=False)

    np.testing.assert_equal(prepared_adata.obsm[basis], test_data.obsm[basis])
    # ensure shape is correct when reading with no keys
    prepared_adata2 = ds_reader.read_dataset(filesystem=fs, path=path, dataset=dataset, keys=dict())
    assert prepared_adata2.shape[0] == test_data.shape[0]


def test_prepare_cxg_tile_db(test_data, measures, dimensions, continuous_obs, basis, tmp_path):
    try:
        from cirrocumulus.tiledb_dataset import TileDBDataset
        output_dir = str(tmp_path)
        test_data = test_data[:, measures]
        test_data.obs = test_data.obs[dimensions + continuous_obs]
        import subprocess
        output_cxg = os.path.join(output_dir, 'test.cxg')
        output_h5ad = os.path.join(output_dir, 'test.h5ad')
        test_data.write(output_h5ad)
        subprocess.check_call(['cellxgene', 'convert', '-o', output_cxg, '--disable-corpora-schema', output_h5ad])
        read_and_diff(TileDBDataset(), output_cxg, test_data, measures, dimensions, continuous_obs, basis)
    except ModuleNotFoundError:
        print('Skipping TileDB tests')


def test_prepare_join_obs_index(test_data, tmp_path):
    output_dir = str(tmp_path)
    test_data2 = test_data.copy()
    test_data2 = test_data2[[0, 1, 2]]
    with pytest.raises(ValueError):
        PrepareData(datasets=[test_data, test_data2], output=output_dir)


@pytest.mark.parametrize("file_format", ['zarr', 'parquet'])
def test_prepare(test_data, measures, dimensions, continuous_obs, basis, file_format, tmp_path):
    file_format2ext = dict(parquet='.cpq', zarr='.zarr')
    output_dir = str(tmp_path / 'test.{}'.format(file_format2ext[file_format]))
    test_data = test_data[:, measures]
    test_data.obs = test_data.obs[dimensions + continuous_obs]
    prepare_data = PrepareData(datasets=[test_data], output=output_dir, output_format=file_format)
    prepare_data.execute()
    if file_format == 'parquet':
        reader = ParquetDataset()
    elif file_format == 'zarr':
        reader = ZarrDataset()
    read_and_diff(reader, output_dir, test_data, measures,
                  dimensions, continuous_obs, basis)


def test_prepare_jsonl(test_data, measures, dimensions, continuous_obs, basis, tmp_path):
    output_dir = str(tmp_path)
    test_data = test_data[:, measures]
    test_data.obs = test_data.obs[dimensions + continuous_obs]
    prepare_data = PrepareData(datasets=[test_data],
                               output=os.path.join(output_dir, 'test.jsonl'), output_format='jsonl')
    prepare_data.execute()
