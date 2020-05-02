import json
import os

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import scipy.sparse

from cirrocumulus.abstract_dataset import AbstractDataset
from cirrocumulus.simple_data import SimpleData


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
            df = pq.read_table(file_system.open(path + '/index.parquet')).to_pandas()

            for column in df:
                result_df[column] = df[column]

        for key in var_keys + obs_keys:
            df = pq.read_table(file_system.open(path + '/' + key + '.parquet')).to_pandas()
            if rename:
                for column in df:
                    result_df['{}_{}'.format(key, column)] = df[column]
            else:
                for column in df:
                    result_df[column] = df[column]
        return result_df

    def read_data_sparse(self, file_system, path, obs_keys=[], var_keys=[], basis=[], dataset=None, schema=None):
        dataset_id = dataset.id
        # if basis, read index to get bins, x, and y
        X = None

        if len(var_keys) > 0:
            data = []
            row = []
            col = []
            X_path = os.path.join(path, 'X')
            shape = schema['shape']
            for i in range(len(var_keys)):
                key = var_keys[i]
                df = pq.read_table(file_system.open(X_path + '/' + key + '.parquet')).to_pandas()
                data.append(df['value'])
                row.append(df['index'])
                col.append(np.repeat(i, len(df)))

            data = np.concatenate(data)
            row = np.concatenate(row)
            col = np.concatenate(col)
            # X = scipy.sparse.csr_matrix((data, (row, col)), shape=(shape[0], len(var_keys)))
            X = scipy.sparse.csc_matrix((data, (row, col)), shape=(shape[0], len(var_keys)))
        obs = None

        for key in obs_keys:
            obs_path = os.path.join(path, 'obs')
            cache_key = str(dataset_id) + '-' + key
            cached_value = self.cached_data.get(cache_key)
            if cached_value is None:
                df = pq.read_table(file_system.open(obs_path + '/' + key + '.parquet'), columns=['value']).to_pandas()
                # ignore index in obs for now
                cached_value = df['value']
            if obs is None:
                obs = pd.DataFrame()
            obs[key] = cached_value
            self.cached_data[cache_key] = cached_value

        if basis is not None and len(basis) > 0:
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
                    df = pq.read_table(file_system.open(basis_path), columns=columns_to_fetch).to_pandas()
                    cached_value = df
                    self.cached_data[cache_key] = cached_value
                if obs is None:
                    obs = pd.DataFrame()
                for c in columns_to_fetch:
                    obs[c] = cached_value[c]
        return SimpleData(X, obs, pd.DataFrame(index=pd.Index(var_keys)))

    def read(self, file_system, path, obs_keys=[], var_keys=[], basis=None, dataset=None, schema=None):
        dataset_id = dataset.id
        if self.cached_dataset_id != dataset_id:
            self.cached_dataset_id = dataset_id
            self.cached_data = {}
        # path is directory
        if not path.endswith('.parquet'):
            return self.read_data_sparse(file_system, path, obs_keys, var_keys, basis, dataset, schema)
        # dense parquet file
        keys = obs_keys + var_keys

        if basis is not None and len(basis) > 0:
            for b in basis:
                keys += b['coordinate_columns']
        if len(keys) == 0:
            return SimpleData(None, pd.DataFrame(), pd.Index([]))
        df = pq.read_table(file_system.open(path), columns=keys).to_pandas()
        X = df[var_keys].values
        if basis is not None and len(basis) > 0:
            for b in basis:
                obs_keys = obs_keys + b['coordinate_columns']
        obs = df[obs_keys]
        return SimpleData(X, obs, pd.DataFrame(index=pd.Index(var_keys)))

    def schema(self, file_system, path):
        if path.endswith('.json') or path.endswith('.json.gz'):  # precomputed dataset
            return super().schema(file_system, path)
        elif path.endswith('.parquet'):
            metadata = None
            with file_system.open(path, 'rb') as s:
                parquet_file = pq.ParquetFile(s)
                schema = parquet_file.schema.to_arrow_schema()
                result = {'version': '1'}
                result['nObs'] = parquet_file.metadata.num_rows
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
                    result['var'] = list(
                        sorted(result['var'],
                            key=lambda x: ('zzzzz' + x.lower()) if x[0].isdigit() else x.lower()))
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
        else:
            return super().schema(file_system, os.path.join(path, 'index.json.gz'))
