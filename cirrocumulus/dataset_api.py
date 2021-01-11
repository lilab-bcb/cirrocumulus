import os

from .file_system_adapter import FileSystemAdapter


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
        self.fs_adapter = FileSystemAdapter()
        self.default_provider = None
        self.cached_schema = None
        self.cached_dataset_id = None

    def get_dataset_provider(self, path):
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
        dataset_id = dataset['id']
        if self.cached_dataset_id == dataset_id:
            return self.cached_schema
        path = dataset['url']
        provider = self.get_dataset_provider(path)
        schema_dict = provider.schema(self.fs_adapter.get_fs(path), path)
        if 'summary' in dataset:
            schema_dict['summary'] = dataset['summary']
        if 'markers' in schema_dict:
            schema_dict['markers_read_only'] = schema_dict.pop('markers')
        self.cached_schema = schema_dict
        self.cached_dataset_id = dataset['id']
        return schema_dict

    def has_precomputed_stats(self, dataset):
        path = dataset['url']
        provider = self.get_dataset_provider(path)
        return provider.has_precomputed_stats(self.fs_adapter.get_fs(path), path, dataset)

    def read_precomputed_stats(self, dataset, obs_keys=[], var_keys=[]):
        path = dataset['url']
        provider = self.get_dataset_provider(path)
        return provider.read_precomputed_stats(self.fs_adapter.get_fs(path), path, obs_keys=obs_keys, var_keys=var_keys)

    def read_precomputed_grouped_stats(self, dataset, obs_keys=[], var_keys=[]):
        path = dataset['url']
        provider = self.get_dataset_provider(path)
        return provider.read_precomputed_grouped_stats(self.fs_adapter.get_fs(path), path, obs_keys=obs_keys,
            var_keys=var_keys)

    def read_precomputed_basis(self, dataset, obs_keys=[], var_keys=[], basis=None):
        path = dataset['url']
        provider = self.get_dataset_provider(path)
        return provider.read_precomputed_basis(self.fs_adapter.get_fs(path), path, obs_keys=obs_keys, var_keys=var_keys,
            basis=basis)

    def read_dataset(self, dataset, keys=[]):
        path = dataset['url']
        provider = self.get_dataset_provider(path)
        return provider.read_dataset(self.fs_adapter.get_fs(path), path, keys=keys, dataset=dataset,
            schema=self.schema(dataset))
