from .file_system import FileSystem


def get_path(dataset, dataset_key):
    path = dataset['url']
    if dataset_key is not None:
        slash_index = path.rfind('/')
        base_url = path[0:slash_index + 1]
        file_name = path[slash_index + 1:]
        file_name = file_name[0:file_name.rfind('.')] + '_' + dataset_key
        path = base_url + file_name + '.parquet'
    return path


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
        if 'summary' in dataset:
            value['summary'] = dataset['summary']
        return value

    def statistics(self, dataset, keys, basis):
        path = dataset['url']
        provider = self.suffix_to_provider[path[path.rfind('.') + 1:].lower()]
        return provider.statistics(self.fs, path, keys, basis)

    def table(self, dataset, keys, basis=None, dataset_key=None):
        path = get_path(dataset, dataset_key)
        provider = self.suffix_to_provider[path[path.rfind('.') + 1:].lower()]
        return provider.table(self.fs, path, keys, basis=basis)

    def tables(self, dataset, keys, basis=None):
        path = dataset['url']
        provider = self.suffix_to_provider[path[path.rfind('.') + 1:].lower()]
        return provider.tables(self.fs, path, keys, basis=basis)
