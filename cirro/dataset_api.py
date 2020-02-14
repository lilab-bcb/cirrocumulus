import os

from .file_system import FileSystem


def get_path(dataset, dataset_path):
    path = dataset['url']
    if path[len(path) - 1] == '/':  # remove trailing slash
        path = path[0:len(path) - 1]
    if path.endswith('.json'):
        path = os.path.dirname(path)
    path = path + '/' + dataset_path
    return path


class DatasetAPI:
    def __init__(self):
        self.suffix_to_provider = {}
        self.fs = FileSystem()
        self.default_provider = None

    def get_provider(self, path):
        index = path.rfind('.')
        if index == -1:
            return self.default_provider
        suffix = path[index + 1:].lower()

        provider = self.suffix_to_provider.get(suffix)
        return provider if provider is not None else self.default_provider

    def add(self, provider):
        if self.default_provider is None:
            self.default_provider = provider
        suffixes = provider.get_suffixes()
        for suffix in suffixes:
            self.suffix_to_provider[suffix.lower()] = provider

    def schema(self, dataset):
        path = dataset['url']
        provider = self.get_provider(path)
        value = provider.schema(self.fs, path)
        if 'summary' in dataset:
            value['summary'] = dataset['summary']
        return value

    def has_precomputed_stats(self, dataset):
        path = dataset['url']
        provider = self.get_provider(path)
        return provider.has_precomputed_stats(self.fs, path, dataset)

    def read_precomputed_stats(self, dataset, obs_keys=[], var_keys=[]):
        path = dataset['url']
        provider = self.get_provider(path)
        return provider.read_precomputed_stats(self.fs, path, obs_keys=obs_keys, var_keys=var_keys)

    def read_precomputed_grouped_stats(self, dataset, obs_keys=[], var_keys=[]):
        path = dataset['url']
        provider = self.get_provider(path)
        return provider.read_precomputed_grouped_stats(self.fs, path, obs_keys=obs_keys, var_keys=var_keys)

    def read_precomputed_basis(self, dataset, obs_keys=[], var_keys=[], basis=None):
        path = dataset['url']
        provider = self.get_provider(path)
        return provider.read_precomputed_basis(self.fs, path, obs_keys=obs_keys, var_keys=var_keys, basis=basis)

    def read_summarized(self, dataset, obs_keys=[], var_keys=[], index=False, rename=False,
                        path=None):
        path = get_path(dataset, path)
        provider = self.get_provider(path)
        return provider.read_summarized(self.fs, path, obs_keys=obs_keys, var_keys=var_keys, index=index,
            rename=rename, dataset=dataset)

    def read(self, dataset, obs_keys=[], var_keys=[], basis=None):
        path = dataset['url']
        provider = self.get_provider(path)
        return provider.read(self.fs, path, obs_keys=obs_keys, var_keys=var_keys, basis=basis, dataset=dataset)
