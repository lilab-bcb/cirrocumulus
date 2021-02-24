import json
import os

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import scipy.sparse

from cirrocumulus.abstract_dataset import AbstractDataset


class ParquetDataset(AbstractDataset):

    def __init__(self):
        super().__init__()
        self.cached_dataset_id = None
        self.cached_data = {}

    def get_suffixes(self):
        return ['parquet', 'pq']

    def read_summarized(self, file_system, path, obs_keys=[], var_keys=[], index=False, rename=False, dataset=None):
        result_df = pd.DataFrame()
        if index:
            df = pq.read_table(path + '/index.parquet', filesystem=file_system).to_pandas()

            for column in df:
                result_df[column] = df[column]

        for key in var_keys + obs_keys:
            df = pq.read_table(path + '/' + key + '.parquet', filesystem=file_system).to_pandas()
            if rename:
                for column in df:
                    result_df['{}_{}'.format(key, column)] = df[column]
            else:
                for column in df:
                    result_df[column] = df[column]
        return result_df

    def read_data_sparse(self, file_system, path, keys, dataset=None, schema=None):
        dataset_id = dataset['id']
        # if basis, read index to get bins, x, and y
        result_df = None
        var_keys = keys.pop('X', [])
        obs_keys = keys.pop('obs', [])
        basis = keys.pop('basis', [])

        if len(var_keys) > 0:
            data = []
            row = []
            col = []
            node_path = os.path.join(path, 'X')
            shape = schema['shape']
            for i in range(len(var_keys)):
                key = var_keys[i]
                df = pq.read_table(node_path + '/' + key + '.parquet', filesystem=file_system).to_pandas()
                data.append(df['value'])
                row.append(df['index'])
                col.append(np.repeat(i, len(df)))

            data = np.concatenate(data)
            row = np.concatenate(row)
            col = np.concatenate(col)
            # X = scipy.sparse.csr_matrix((data, (row, col)), shape=(shape[0], len(var_keys)))
            X = scipy.sparse.csc_matrix((data, (row, col)), shape=(shape[0], len(var_keys)))
            result_df = pd.DataFrame.sparse.from_spmatrix(X, columns=var_keys)

        for node in keys.keys():
            if result_df is None:
                result_df = pd.DataFrame()
            features = keys[node]
            node_path = os.path.join(path, node)
            for i in range(len(features)):
                key = features[i]
                df = pq.read_table(node_path + '/' + key + '.parquet', filesystem=file_system).to_pandas()
                result_df[key] = df['value']
        if len(obs_keys) > 0:
            if result_df is None:
                result_df = pd.DataFrame()
            for key in obs_keys:
                node_path = os.path.join(path, 'obs')
                cache_key = str(dataset_id) + '-' + key
                cached_value = self.cached_data.get(cache_key)
                if cached_value is None:
                    df = pq.read_table(node_path + '/' + key + '.parquet', filesystem=file_system,
                        columns=['value']).to_pandas()
                    # ignore index in obs for now
                    cached_value = df['value']
                result_df[key] = cached_value
                self.cached_data[cache_key] = cached_value

        if basis is not None and len(basis) > 0:
            if result_df is None:
                result_df = pd.DataFrame()
            obsm_path = os.path.join(path, 'obsm')
            for b in basis:
                cache_key = str(dataset_id) + '-' + b['full_name']
                is_precomputed = b['precomputed']
                if is_precomputed:
                    columns_to_fetch = [b['full_name']]
                else:  # need coordinates
                    columns_to_fetch = b['coordinate_columns']
                cached_value = self.cached_data.get(cache_key)
                if cached_value is None:
                    basis_path = (obsm_path + '/') + (b['full_name' if is_precomputed else 'name']) + '.parquet'
                    df = pq.read_table(basis_path, filesystem=file_system, columns=columns_to_fetch).to_pandas()
                    cached_value = df
                    self.cached_data[cache_key] = cached_value

                for c in columns_to_fetch:
                    result_df[c] = cached_value[c]
        if result_df is None:
            shape = schema['shape']
            result_df = pd.DataFrame(index=pd.RangeIndex(shape[0]))
        return result_df

    def read_data_dense(self, file_system, path, keys=None, dataset=None, schema=None):
        # dense parquet file
        var_keys = keys.pop('X', [])
        obs_keys = keys.pop('obs', [])
        basis = keys.pop('basis', [])
        all_keys = obs_keys + var_keys
        for key in keys.keys():
            all_keys += keys[key]

        if basis is not None and len(basis) > 0:
            for b in basis:
                all_keys += b['coordinate_columns']
        if len(all_keys) == 0:
            return pq.read_table(path, filesystem=file_system, columns=['index']).to_pandas()
        return pq.read_table(path, filesystem=file_system, columns=all_keys).to_pandas()

    def read_dataset(self, file_system, path, keys=None, dataset=None, schema=None):
        dataset_id = dataset['id']
        if self.cached_dataset_id != dataset_id:
            self.cached_dataset_id = dataset_id
            self.cached_data = {}
        if keys is None:
            keys = {}
        # path is directory
        keys = keys.copy()
        if not path.endswith('.parquet'):
            return self.read_data_sparse(file_system, path, keys, dataset, schema)

        return self.read_data_dense(file_system, path, keys, dataset, schema)

    def schema(self, file_system, path):
        if path.endswith('.json') or path.endswith('.json.gz'):
            return super().schema(file_system, path)
        elif path.endswith('.parquet'):
            metadata = None
            with file_system.open(path, 'rb') as s:
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
            return super().schema(file_system, os.path.join(path, 'index.json.gz'))
