import json
from abc import abstractmethod, ABC


class AbstractDataset(ABC):

    def __init__(self):
        super().__init__()
        self.cached_dataset_id = None
        self.cached_data = {}

    @abstractmethod
    def get_suffixes(self):
        pass

    def has_precomputed_stats(self, file_system, path, dataset):
        return False

    @abstractmethod
    def read(self, file_system, path, obs_keys=[], var_keys=[], basis=None, dataset=None):
        pass

    def schema(self, file_system, path):
        if path.endswith('.gz'):
            import gzip
            with gzip.open(file_system.open(path)) as s:
                return json.load(s)
        else:
            with file_system.open(path) as s:
                return json.load(s)
