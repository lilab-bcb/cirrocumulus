import json
import os

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import scipy.sparse
import sys

from cirro.simple_data import SimpleData


class ParquetDataset:

    def __init__(self):
        self.cached_path = None
        self.cached_parquet_file = None
        self.cached_stream = None
        self.cached_dataset_id = None
        self.cached_data = {}

    def get_suffixes(self):
        return ['parquet', 'pq', 'json']

    def get_cached_stream(self, file_system, path):
        if self.cached_path != path:
            self.cached_path = path
            if self.cached_stream is not None:
                self.cached_stream.close()
            self.cached_stream = file_system.open(path)
            return self.cached_stream
        return None

    def get_file(self, file_system, path):
        stream = self.get_cached_stream(file_system, path)
        if stream is not None:  # stream updated
            self.cached_parquet_file = pq.ParquetFile(self.cached_stream)
        return self.cached_parquet_file

    def close(self, path):
        if self.cached_path == path and self.cached_stream is not None:
            self.cached_stream.close()

    def schema(self, file_system, path):
        if path.endswith('.json'):  # prepared dataset
            with file_system.open(path) as s:
                return json.load(s)
        self.get_cached_stream(file_system, path)
        parquet_file = self.get_file(file_system, path)
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

    def statistics(self, file_system, path, keys, basis):
        keys = ParquetDataset.get_keys(keys, basis=basis)
        parquet_file = self.get_file(file_system, path)
        schema = parquet_file.schema.to_arrow_schema()
        indices = []
        key_to_stats = {}

        for key in keys:
            indices.append(schema.get_field_index(key))
            min_max = [sys.float_info.max, sys.float_info.min]
            key_to_stats[key] = min_max

        for i in range(parquet_file.num_row_groups):
            row_group = parquet_file.metadata.row_group(i)
            for j in range(len(indices)):
                min_max = key_to_stats[keys[j]]
                stats = row_group.column(indices[j]).statistics
                min_max[0] = min(stats.min, min_max[0])
                min_max[1] = max(stats.max, min_max[1])
        return key_to_stats

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

    def read_data_sparse(self, file_system, path, obs_keys=[], var_keys=[], basis=None, dataset=None):
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
                df = pq.read_table(file_system.open(data_path + '/' + key + '.parquet')).to_pandas()
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
                df = pq.read_table(file_system.open(data_path + '/' + key + '.parquet'),
                    columns=['value']).to_pandas()  # ignore index in obs for now
                cached_value = df['value']
            if obs is None:
                obs = pd.DataFrame()
            obs[key] = cached_value
            self.cached_data[cache_key] = cached_value

        if basis is not None:
            # don't need the coordinates, only bins
            cache_key = str(dataset_id) + '-' + basis['name']
            cached_value = self.cached_data.get(cache_key)
            if cached_value is None:
                df = pq.read_table(file_system.open(data_path + '/' + basis['name'] + '.parquet'),
                    columns=['index']).to_pandas()
                cached_value = df['index']
            self.cached_data[cache_key] = cached_value
            if obs is None:
                obs = pd.DataFrame(index=cached_value)
            else:
                obs.index = cached_value

        return SimpleData(X, obs, pd.DataFrame(index=pd.Index(var_keys)))

    def read(self, file_system, path, obs_keys=[], var_keys=[], basis=None, dataset=None):
        # path is path to index.json
        if path.endswith('.json'):
            return self.read_data_sparse(file_system, path, obs_keys, var_keys, basis, dataset)
        parquet_file = self.get_file(file_system, path)
        keys = ParquetDataset.get_keys(obs_keys + var_keys, basis=basis)
        if len(keys) == 0:
            return SimpleData(None, pd.DataFrame(index=pd.RangeIndex(parquet_file.metadata.num_rows)),
                pd.Index([]), {})
        df = pq.read_table(parquet_file, keys).to_pandas()
        X = df[var_keys].values
        obs = df[obs_keys]

        if basis is not None:
            for key in basis['coordinate_columns']:
                df[obs] = df[key]
        return SimpleData(X, obs, pd.DataFrame(index=pd.Index(var_keys)))
