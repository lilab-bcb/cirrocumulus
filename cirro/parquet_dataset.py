import json

import pyarrow as pa
import pyarrow.parquet as pq

from cirro.abstract_dataset import AbstractDataset


class ParquetDataset(AbstractDataset):

    def __init__(self):
        super().__init__()

    def get_suffixes(self):
        return ['parquet', 'pq', 'json', 'pjson']

    def to_pandas(self, file_system, path, columns=None):
        return pq.read_table(file_system.open(path + '.parquet'), columns=columns).to_pandas()

    def schema(self, file_system, path):
        if path.endswith('.json'):  # precomputed dataset
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
