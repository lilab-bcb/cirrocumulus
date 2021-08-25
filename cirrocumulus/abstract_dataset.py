import json
from abc import abstractmethod, ABC


class AbstractDataset(ABC):

    def __init__(self):
        super().__init__()

    @abstractmethod
    def get_suffixes(self):
        pass

    @abstractmethod
    def read_dataset(self, filesystem, path, keys=None, dataset=None):
        pass

    def schema(self, filesystem, path):
        if path.endswith('.gz'):
            import gzip
            with gzip.open(filesystem.open(path)) as s:
                return json.load(s)
        else:
            with filesystem.open(path) as s:
                return json.load(s)
