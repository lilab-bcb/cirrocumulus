import os

import json
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import scipy.sparse

from cirro.simple_data import SimpleData


class ParquetDataset:

    def __init__(self):
        self.cached_dataset_id = None
        self.cached_data = {}
        self.cached_schema = None
        self.cached_schema_path = None

    def get_suffixes(self):
        return ['parquet', 'pq', 'json']

    def schema(self, file_system, path):
        if path == self.cached_schema_path:
            return self.cached_schema
        if path.endswith('.json'):  # precomputed dataset
            with file_system.open(path) as s:
                s = json.load(s)
                self.cached_schema = s
                self.cached_schema_path = path
                return s
        with file_system.open(path, 'rb') as s:
            parquet_file = pq.ParquetFile(s)
            schema = parquet_file.schema.to_arrow_schema()
            metadata = json.loads(schema.metadata[b'pegasus'])
            result = {'version': '1'}
            all_obs = metadata['obs']
            obs = []
            obs_cat = []
            for name in all_obs:
                field = schema.field(name)
                if isinstance(field.type, pa.lib.DictionaryType):
                    obs_cat.append(name)
                else:
                    obs.append(name)
            result['var'] = metadata['var']
            result['obs'] = obs
            result['obsCat'] = obs_cat
            result['nObs'] = parquet_file.metadata.num_rows
            result['embeddings'] = metadata['obsm']
            return result

    @staticmethod
    def get_keys(keys, basis=None):
        if basis is not None:
            keys = keys + basis['coordinate_columns']
        return keys

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

    def read_data_sparse(self, file_system, path, obs_keys=[], var_keys=[], basis=[], dataset=None):
        # path is path to index.json
        # if path ends with /data then X is stored as index, value pairs
        # if basis, read index to get bins, x, and y
        path = os.path.dirname(path)
        data_path = os.path.join(path, 'data')
        data = []
        row = []
        col = []
        X = None
        nobs = dataset['nObs']

        if len(var_keys) > 0:
            for i in range(len(var_keys)):
                key = var_keys[i]
                with file_system.open(data_path + '/' + key + '.parquet', 'rb') as s:
                    df = pq.read_table(s).to_pandas()
                data.append(df['value'])
                row.append(df['index'])
                col.append(np.repeat(i, len(df)))

            data = np.concatenate(data)
            row = np.concatenate(row)
            col = np.concatenate(col)
            X = scipy.sparse.csr_matrix((data, (row, col)), shape=(nobs, len(var_keys)))
        obs = None
        dataset_id = dataset.id
        if self.cached_dataset_id != dataset_id:
            self.cached_dataset_id = dataset_id
            self.cached_data = {}

        for key in obs_keys:
            cache_key = str(dataset_id) + '-' + key
            cached_value = self.cached_data.get(cache_key)
            if cached_value is None:
                with file_system.open(data_path + '/' + key + '.parquet', 'rb') as s:
                    df = pq.read_table(s, columns=['value']).to_pandas()  # ignore index in obs for now
                cached_value = df['value']
            if obs is None:
                obs = pd.DataFrame()
            obs[key] = cached_value
            self.cached_data[cache_key] = cached_value

        if basis is not None and len(basis) > 0:
            for b in basis:
                cache_key = str(dataset_id) + '-' + b['full_name']
                columns_to_fetch = [b['full_name']]
                if not b['precomputed']:  # need coordinates and bins
                    columns_to_fetch += b['coordinate_columns']
                cached_value = self.cached_data.get(cache_key)
                if cached_value is None:
                    with file_system.open(data_path + '/' + b['full_name'] + '.parquet', 'rb') as s:
                        df = pq.read_table(s, columns=columns_to_fetch).to_pandas()
                    cached_value = df
                    self.cached_data[cache_key] = cached_value
                if obs is None:
                    obs = pd.DataFrame()
                for c in columns_to_fetch:
                    obs[c] = cached_value[c]
        return SimpleData(X, obs, pd.DataFrame(index=pd.Index(var_keys)))

    def read(self, file_system, path, obs_keys=[], var_keys=[], basis=None, dataset=None):
        # path is path to index.json
        if path.endswith('.json'):
            return self.read_data_sparse(file_system, path, obs_keys, var_keys, basis, dataset)

        keys = ParquetDataset.get_keys(obs_keys + var_keys, basis=basis)
        if len(keys) == 0:
            return SimpleData(None, pd.DataFrame(), pd.Index([]))
        with file_system.open(path, 'rb') as s:
            df = pq.read_table(s, keys).to_pandas()
        X = df[var_keys].values
        if basis is not None and len(basis) > 0:
            for b in basis:
                obs_keys = obs_keys + b['coordinate_columns']
        obs = df[obs_keys]
        return SimpleData(X, obs, pd.DataFrame(index=pd.Index(var_keys)))
