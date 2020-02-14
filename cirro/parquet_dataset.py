import json
import os

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from cirro.abstract_dataset import AbstractDataset
from cirro.simple_data import SimpleData


def write_pq(d, output_dir, name, write_statistics=True, row_group_size=None):
    os.makedirs(output_dir, exist_ok=True)
    pq.write_table(pa.Table.from_pydict(d), output_dir + os.path.sep + name + '.parquet',
        write_statistics=write_statistics, row_group_size=row_group_size)


class ParquetDataset(AbstractDataset):

    def __init__(self):
        super().__init__()

    def get_suffixes(self):
        return ['parquet', 'pq', 'json', 'pjson']


    def read_data_sparse(self, file_system, path, obs_keys=[], var_keys=[], basis=[], dataset=None):
        # path is path to index.json
        # if path ends with /data then X is stored as index, value pairs
        # if basis, read index to get bins, x, and y
        path = os.path.dirname(path)
        data_path = os.path.join(path, 'data')
        dataset_attrs = self.get_dataset_attrs(file_system, path)
        shape = dataset_attrs['shape']

        data = []
        row = []
        col = []
        X = None

        if len(var_keys) > 0:
            for i in range(len(var_keys)):
                key = var_keys[i]
                df = self.to_pandas(file_system, data_path + '/' + key)
                data.append(df['value'])
                row.append(df['index'])
                col.append(np.repeat(i, len(df)))

            data = np.concatenate(data)
            row = np.concatenate(row)
            col = np.concatenate(col)
            X = scipy.sparse.csr_matrix((data, (row, col)), shape=(shape[0], len(var_keys)))
        obs = None
        dataset_id = dataset.id
        if self.cached_dataset_id != dataset_id:
            self.cached_dataset_id = dataset_id
            self.cached_data = {}

        for key in obs_keys:
            cache_key = str(dataset_id) + '-' + key
            cached_value = self.cached_data.get(cache_key)
            if cached_value is None:
                df = self.to_pandas(file_system, data_path + '/' + key, ['value'])  # ignore index in obs for now
                cached_value = df['value']
            if obs is None:
                obs = pd.DataFrame()
            obs[key] = cached_value
            self.cached_data[cache_key] = cached_value

        if basis is not None and len(basis) > 0:
            for b in basis:
                cache_key = str(dataset_id) + '-' + b['full_name']
                is_precomputed = b['precomputed']
                if is_precomputed:
                    columns_to_fetch = [b['full_name']]
                else:  # need coordinates
                    columns_to_fetch = b['coordinate_columns']
                cached_value = self.cached_data.get(cache_key)
                if cached_value is None:
                    basis_path = (data_path + '/') + (b['full_name' if is_precomputed else 'name'])
                    df = self.to_pandas(file_system, basis_path, columns_to_fetch)
                    cached_value = df
                    self.cached_data[cache_key] = cached_value
                if obs is None:
                    obs = pd.DataFrame()
                for c in columns_to_fetch:
                    obs[c] = cached_value[c]
        return SimpleData(X, obs, pd.DataFrame(index=pd.Index(var_keys)))

    @staticmethod
    def get_keys(keys, basis=None):
        if basis is not None:
            keys = keys + basis['coordinate_columns']
        return keys

    def read_summarized(self, file_system, path, obs_keys=[], var_keys=[], index=False, rename=False, dataset=None):
        result_df = pd.DataFrame()

        if index:
            df = self.to_pandas(file_system, path + '/index')
            for column in df:
                result_df[column] = df[column]

        for key in var_keys + obs_keys:
            df = self.to_pandas(file_system, path + '/' + key)
            if rename:
                for column in df:
                    result_df['{}_{}'.format(key, column)] = df[column]
            else:
                for column in df:
                    result_df[column] = df[column]
        return result_df

    def has_precomputed_stats(self, file_system, path, dataset):
        return path.endswith('json')

    def read(self, file_system, path, obs_keys=[], var_keys=[], basis=None, dataset=None):
        # path is path to index.json
        if path.endswith('json'):
            return self.read_data_sparse(file_system, path, obs_keys, var_keys, basis, dataset)
        keys = []
        if basis is not None and len(basis) > 0:
            for b in basis:
                keys += AbstractDataset.get_keys(obs_keys + var_keys, basis=b)
        if len(keys) == 0:
            return SimpleData(None, pd.DataFrame(), pd.Index([]))
        df = self.to_pandas(file_system, path, keys)
        X = df[var_keys].values
        if basis is not None and len(basis) > 0:
            for b in basis:
                obs_keys = obs_keys + b['coordinate_columns']
        obs = df[obs_keys]
        return SimpleData(X, obs, pd.DataFrame(index=pd.Index(var_keys)))


    def to_pandas(self, file_system, path, columns=None):
        return pq.read_table(file_system.open(path + '.parquet'), columns=columns).to_pandas()

    def schema(self, file_system, path):
        if path.endswith('.json') or path.endswith('.json.gz'):  # precomputed dataset
            return super().schema(file_system, path)
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
