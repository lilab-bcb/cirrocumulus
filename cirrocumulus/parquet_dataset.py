import concurrent.futures
import json
import os

import anndata
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import scipy.sparse
from anndata import AnnData

from cirrocumulus.abstract_dataset import AbstractDataset

max_workers = min(12, pa.cpu_count())
executor = concurrent.futures.ThreadPoolExecutor(max_workers=max_workers)


def read_table(path, filesystem, columns=None):
    return pq.read_table(path, filesystem=filesystem, columns=columns, use_threads=False)


def read_df(path, filesystem, columns=None):
    return read_table(path, filesystem=filesystem, columns=columns)


def read_tables(paths, filesystem, columns=None):
    futures = []
    for path in paths:
        future = executor.submit(read_table, path, filesystem=filesystem, columns=columns)
        futures.append(future)
    concurrent.futures.wait(futures)
    return futures


def read_dfs(paths, filesystem, columns=None):
    futures = []
    for path in paths:
        future = executor.submit(read_df, path, filesystem=filesystem, columns=columns)
        futures.append(future)
    concurrent.futures.wait(futures)
    return futures


def get_matrix(futures, shape=None):
    data = []
    row = []
    col = []
    is_sparse = None
    for i in range(len(futures)):
        t = futures[i].result()
        if i == 0:
            is_sparse = 'index' in t.column_names
        if is_sparse:
            row.append(t.column('index').to_numpy())
            col.append(np.repeat(i, len(t)))
        data.append(t.column('value').to_numpy())
    if len(data) == 0:
        raise ValueError('No data')
    if len(row) > 0:  # sparse
        data = np.concatenate(data)
        row = np.concatenate(row)
        col = np.concatenate(col)
        # X = scipy.sparse.coo_matrix((data, (row, col)), shape=(shape[0], len(keys))).to_csc()
        # X = scipy.sparse.csr_matrix((data, (row, col)), shape=(shape[0], len(var_keys)))
        X = scipy.sparse.csc_matrix((data, (row, col)), shape=(shape[0], len(futures)))
    else:
        X = np.array(data).T
    return X


class ParquetDataset(AbstractDataset):

    def __init__(self):
        super().__init__()

    def get_suffixes(self):
        return ['parquet', 'pq', 'cpq']

    def read_data_sparse(self, filesystem, path, keys, dataset=None, schema=None):
        shape = schema['shape']
        schema_version = schema.get('version', '1.0.0')

        X = None
        adata_modules = None
        obs = None
        var = None
        obsm = {}
        var_keys = keys.pop('X', [])
        obs_keys = keys.pop('obs', [])
        basis = keys.pop('basis', [])
        module_keys = keys.pop('module', [])
        if len(var_keys) > 0:
            node_path = os.path.join(path, 'X')
            paths = [node_path + '/' + key + '.parquet' for key in var_keys]
            X = get_matrix(read_tables(paths, filesystem), shape)
            var = pd.DataFrame(index=var_keys)
        if len(module_keys) > 0:
            node_path = os.path.join(path, 'X_module')
            paths = [node_path + '/' + key + '.parquet' for key in module_keys]
            module_X = get_matrix(read_tables(paths, filesystem))
            adata_modules = anndata.AnnData(X=module_X, var=pd.DataFrame(index=module_keys))
        if len(obs_keys) > 0:
            obs = pd.DataFrame()
            node_path = os.path.join(path, 'obs')
            paths = [node_path + '/' + key + '.parquet' for key in obs_keys]
            futures = read_tables(paths, filesystem, columns=['value'])
            for i in range(len(futures)):
                df = futures[i].result().to_pandas()
                obs[obs_keys[i]] = df['value']

        if basis is not None and len(basis) > 0:
            node_path = os.path.join(path, 'obsm')
            paths = [node_path + '/' + b['name'] + '.parquet' for b in basis]
            futures = read_tables(paths, filesystem)
            for i in range(len(futures)):
                table = futures[i].result()
                b = basis[i]
                vals = []
                columns_to_fetch = b['coordinate_columns']
                for c in columns_to_fetch:
                    vals.append(table.column(c))
                obsm[b['name']] = np.array(vals).T
        adata = AnnData(X=X, obs=obs, var=var, obsm=obsm)
        if adata_modules is not None:
            adata.uns['X_module'] = adata_modules
        return adata

    def read_data_dense(self, filesystem, path, keys=None, dataset=None, schema=None):
        var_keys = keys.pop('X', [])
        obs_keys = keys.pop('obs', [])
        basis = keys.pop('basis', [])
        all_keys = obs_keys + var_keys
        X = None
        obs = None
        var = None
        obsm = {}
        for key in keys.keys():
            all_keys += keys[key]

        if len(basis) > 0:
            for b in basis:
                all_keys += b['coordinate_columns']
        df = read_table(path, filesystem=filesystem, columns=all_keys)
        if len(var_keys) > 0:
            X = df[var_keys]
        if len(obs_keys) > 0:
            obs = df[obs_keys]
        if len(basis) > 0:
            for b in basis:
                obsm[b['name']] = df[b['coordinate_columns']]
        return AnnData(X=X, obs=obs, var=var, obsm=obsm)

    def read_dataset(self, filesystem, path, keys=None, dataset=None, schema=None):
        if keys is None:
            keys = {}
        # path is directory
        keys = keys.copy()
        if not path.endswith('.parquet'):
            return self.read_data_sparse(filesystem, path, keys, dataset, schema)
        return self.read_data_dense(filesystem, path, keys, dataset, schema)

    def schema(self, filesystem, path):
        if path.endswith('.json') or path.endswith('.json.gz'):
            return super().schema(filesystem, path)
        elif path.endswith('.parquet'):
            metadata = None
            with filesystem.open(path, 'rb') as s:
                parquet_file = pq.ParquetFile(s)
                schema = parquet_file.schema.to_arrow_schema()
                result = {'version': '1'}
                for metadata_key in [b'pegasus']:
                    if metadata_key in schema.metadata:
                        metadata = json.loads(schema.metadata[metadata_key])
                        break
                obs = []
                obs_cat = []
                if metadata is None:
                    # obs, the obsm, then var
                    names = parquet_file.schema.names
                    section = 'obs'
                    result['embeddings'] = []
                    result['var'] = []
                    embedding_basename_to_count = {}
                    for name in names:
                        if section == 'obs' and name.startswith('X_'):
                            section = 'obsm'
                        if section == 'obs':
                            field = schema.field(name)
                            if isinstance(field.type, pa.lib.DictionaryType):
                                obs_cat.append(name)
                            else:
                                obs.append(name)
                        if section == 'obsm' and not name.startswith('X_'):
                            section = 'var'
                        if section == 'obsm':  # X_fitsne_1', 'X_fitsne_2, etc
                            basename = name[0:name.rindex('_')]
                            count = embedding_basename_to_count.get(basename, 0)
                            embedding_basename_to_count[basename] = count + 1
                        if section == 'var':
                            result['var'].append(name)
                    for name in embedding_basename_to_count.keys():
                        result['embeddings'].append(dict(name=name, dimensions=2))
                        if embedding_basename_to_count[name] == 3:
                            result['embeddings'].append(dict(name=name, dimensions=3))
                else:
                    all_obs = metadata['obs']

                    for name in all_obs:
                        field = schema.field(name)
                        if isinstance(field.type, pa.lib.DictionaryType):
                            obs_cat.append(name)
                        else:
                            obs.append(name)
                    result['var'] = metadata['var']
                    result['embeddings'] = metadata['obsm']

                result['obs'] = obs
                result['obsCat'] = obs_cat
                result['shape'] = (parquet_file.metadata.num_rows, len(result['var']))
                return result
        else:  # directory
            return super().schema(filesystem, os.path.join(path, 'index.json.gz'))
