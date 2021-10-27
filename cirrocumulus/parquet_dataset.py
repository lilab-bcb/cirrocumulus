import concurrent.futures
import json
import os

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


def read_tables(paths, filesystem, columns=None):
    futures = []
    for path in paths:
        future = executor.submit(read_table, path, filesystem=filesystem, columns=columns)
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

    def read_data_sparse(self, filesystem, path, keys, dataset=None):
        dataset_info = self.get_dataset_info(filesystem, path)
        shape = dataset_info['shape']
        X = None
        obs = None
        var = None
        obsm = {}
        X_keys = keys.pop('X', [])
        obs_keys = keys.pop('obs', [])
        basis_keys = keys.pop('basis', [])

        if len(X_keys) > 0:
            if len(X_keys) == 1 and isinstance(X_keys[0], slice):  # special case if slice specified for performance
                get_item_x = X_keys[0]
                X_keys = dataset_info['var'][get_item_x]
            node_path = os.path.join(path, 'X')
            paths = [node_path + '/' + key + '.parquet' for key in X_keys]
            X = get_matrix(read_tables(paths, filesystem), shape)
            var = pd.DataFrame(index=X_keys)
        if len(obs_keys) > 0:
            obs = pd.DataFrame()
            node_path = os.path.join(path, 'obs')
            paths = [node_path + '/' + key + '.parquet' for key in obs_keys]
            futures = read_tables(paths, filesystem, columns=['value'])
            for i in range(len(futures)):
                df = futures[i].result().to_pandas()
                obs[obs_keys[i]] = df['value']

        if len(basis_keys) > 0:
            node_path = os.path.join(path, 'obsm')
            paths = [node_path + '/' + key + '.parquet' for key in basis_keys]
            futures = read_tables(paths, filesystem)
            for i in range(len(futures)):
                table = futures[i].result()
                vals = []
                for c in table.column_names:
                    vals.append(table.column(c))
                vals = np.array(vals).T
                obsm[basis_keys[i]] = vals
                if X is None:
                    X = scipy.sparse.coo_matrix((vals.shape[0], 0))
        if X is None and obs is None and len(obsm.keys()) == 0:
            obs = pd.DataFrame(index=pd.RangeIndex(shape[0]).astype(str))
        if X is None and obs is None and len(obsm.keys()) == 0:
            X = scipy.sparse.csc_matrix(tuple(shape))
        adata = AnnData(X=X, obs=obs, var=var, obsm=obsm)
        return adata

    def read_data_dense(self, filesystem, path, keys=None, dataset=None):
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
            for key in basis:
                # 2d only
                all_keys += ['{}_{}'.format(key, 1), '{}_{}'.format(key, 2)]
        df = read_table(path, filesystem=filesystem, columns=all_keys)
        if len(var_keys) > 0:
            X = df[var_keys]
        if len(obs_keys) > 0:
            obs = df[obs_keys]
        if len(basis) > 0:
            for key in basis:
                obsm[key] = df[['{}_{}'.format(key, 1), '{}_{}'.format(key, 2)]]
        return AnnData(X=X, obs=obs, var=var, obsm=obsm)

    def read_dataset(self, filesystem, path, keys=None, dataset=None):
        if keys is None:
            keys = {}
        # path is directory
        keys = keys.copy()
        if not path.endswith('.parquet'):
            return self.read_data_sparse(filesystem, path, keys, dataset)
        return self.read_data_dense(filesystem, path, keys, dataset)

    def get_schema(self, filesystem, path):
        if path.endswith('.json') or path.endswith('.json.gz'):
            return super().get_schema(filesystem, path)
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
            return super().get_schema(filesystem, os.path.join(path, 'index.json.gz'))
