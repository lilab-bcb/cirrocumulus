import json
import sys

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


class ParquetDataset:

    def __init__(self):
        self.cached_path = None
        self.cached_parquet_file = None
        self.cached_stream = None

    def get_file(self, file_system, path):
        if self.cached_path != path:
            if self.cached_stream is not None:
                self.cached_stream.close()
            self.cached_stream = file_system.open(path)
            self.cached_parquet_file = pq.ParquetFile(self.cached_stream)
            self.cached_path = path
        return self.cached_parquet_file

    def close(self, path):
        if self.cached_path == path:
            if self.cached_stream is not None:
                self.cached_stream.close()

    def schema(self, file_system, path):
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

        obsm = metadata['obsm']
        result['var'] = metadata['var']
        result['obs'] = obs
        result['obs_cat'] = obs_cat
        result['n_obs'] = parquet_file.metadata.num_rows
        result['embeddings'] = obsm
        return result

    @staticmethod
    def get_keys(keys, basis=None, index=False):
        if basis is not None:
            keys = keys + basis['coordinate_columns']
        if index:
            keys = keys + ['index']  # get pandas index
        return keys

    def statistics(self, file_system, path, keys):
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

    def tables(self, file_system, path, keys, basis=None, index=False):
        keys = ParquetDataset.get_keys(keys, basis=basis, index=index)
        parquet_file = self.get_file(file_system, path)

        for i in range(parquet_file.num_row_groups):
            yield parquet_file.read_row_group(i, columns=keys)

    def get_df(self, file_system, path, keys, basis=None, index=False):
        keys = ParquetDataset.get_keys(keys, basis=basis, index=index)
        parquet_file = self.get_file(file_system, path)
        if len(keys) == 0:
            return pd.DataFrame(index=pd.RangeIndex(parquet_file.metadata.num_rows))
        table = parquet_file.read(keys)
        return table.to_pandas()
