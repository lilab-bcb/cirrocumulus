import gzip
import json
import os
from abc import abstractmethod, ABC

import cirrocumulus.diff_exp


class AbstractDataset(ABC):

    def __init__(self):
        super().__init__()

    @abstractmethod
    def get_suffixes(self):
        pass

    @abstractmethod
    def read(self, file_system, path, obs_keys=[], var_keys=[], basis=None, dataset=None):
        pass

    def diff_exp(self, file_system, path, mask, dataset, schema, var_range):
        var_names = schema['var']
        start = var_range[0]
        end = var_range[1]
        if end - start > 5000:
            raise ValueError()
        var_names = var_names[start:end]
        adata = self.read(file_system, path, obs_keys=[], var_keys=var_names, basis=[],
            dataset=dataset, schema=schema)
        result = cirrocumulus.diff_exp.diff_exp(adata.X, mask)
        result['name'] = var_names
        # if result is None:
        #     result = batch_result
        # else:
        #     for key in batch_result:
        #         result[key] = np.concatenate((result[key], batch_result[key]))
        return result

    def has_precomputed_stats(self, file_system, path, dataset):
        return file_system.exists(os.path.join(path, 'stats'))

    def read_precomputed_stats(self, file_system, path, obs_keys=[], var_keys=[]):
        result = {}
        if len(obs_keys) > 0:
            for key in obs_keys:
                with gzip.open(file_system.open(path + '/stats/obs/' + key + '.json.gz', 'rb')) as f:
                    result[key] = json.load(f)

        if len(var_keys) > 0:
            for key in var_keys:
                with gzip.open(file_system.open(path + '/stats/X/' + key + '.json.gz', 'rb')) as f:
                    result[key] = json.load(f)
        return result

    def read_precomputed_grouped_stats(self, file_system, path, obs_keys=[], var_keys=[]):
        results = []
        if len(obs_keys) > 0 and len(var_keys) > 0:

            for obs_key in obs_keys:
                values = []
                with gzip.open(file_system.open(path + '/grouped_stats/obs/' + obs_key + '/index.json.gz', 'rb')) as f:
                    categories = json.load(f)
                result = {'categories': categories, 'name': obs_key, 'values': values}
                results.append(result)
                for var_key in var_keys:
                    with gzip.open(
                            file_system.open(path + '/grouped_stats/obs/' + obs_key + '/X/' + var_key + '.json.gz',
                                'rb')) as f:
                        r = json.load(f)
                        r['name'] = var_key
                        values.append(r)
        return results

    def schema(self, file_system, path):
        if path.endswith('.gz'):
            import gzip
            with gzip.open(file_system.open(path)) as s:
                return json.load(s)
        else:
            with file_system.open(path) as s:
                return json.load(s)
