from .file_system import FileSystem


class DatasetAPI:
    def __init__(self):
        self.suffix_to_provider = {}
        self.fs = FileSystem()

    def add(self, suffixes, provider):
        for suffix in suffixes:
            self.suffix_to_provider[suffix.lower()] = provider

    def schema(self, dataset):
        path = dataset['url']
        provider = self.suffix_to_provider[path[path.rfind('.') + 1:].lower()]
        value = provider.schema(self.fs, path)
        return value

    def statistics(self, dataset, keys):
        path = dataset['url']
        provider = self.suffix_to_provider[path[path.rfind('.') + 1:].lower()]
        return provider.statistics(self.fs, path, keys)

    def tables(self, dataset, keys, basis=None):
        path = dataset['url']
        provider = self.suffix_to_provider[path[path.rfind('.') + 1:].lower()]
        return provider.tables(self.fs, path, keys, basis=basis)
