import json
import os

import pyarrow as pa
import pyarrow.parquet as pq

from cirro.abstract_dataset import AbstractDataset


def write_pq(d, output_dir, name, write_statistics=True, row_group_size=None):
    os.makedirs(output_dir, exist_ok=True)
    pq.write_table(pa.Table.from_pydict(d), output_dir + os.path.sep + name + '.parquet',
        write_statistics=write_statistics, row_group_size=row_group_size)


class ParquetDataset(AbstractDataset):

    def __init__(self):
        super().__init__()

    def get_suffixes(self):
        return ['parquet', 'pq', 'json', 'pjson']

    def to_pandas(self, file_system, path, columns=None):
        return pq.read_table(file_system.open(path + '.parquet'), columns=columns).to_pandas()

    def schema(self, file_system, path):
        if path.endswith('.json') or path.endswith('.json.gz'):  # precomputed dataset
            super().schema(file_system, path)
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
