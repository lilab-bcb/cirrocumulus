import os

from cirrocumulus.util import get_fs


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
        self.default_provider = None
        self.cached_dataset_info = None
        self.cached_dataset_id = None

    def get_dataset_provider(self, path):
        index = path.rfind('.')
        if index == -1:
            return self.default_provider
        suffix = path[index + 1:].lower()
        suffix = suffix.rstrip('/')
        provider = self.suffix_to_provider.get(suffix)
        return provider if provider is not None else self.default_provider

    def add(self, provider):
        if self.default_provider is None:
            self.default_provider = provider
        suffixes = provider.get_suffixes()
        for suffix in suffixes:
            self.suffix_to_provider[suffix.lower()] = provider

    def get_dataset_info(self, dataset):
        dataset_id = dataset['id']
        if self.cached_dataset_id == dataset_id:
            dataset_info = self.cached_dataset_info
        else:
            path = dataset['url']
            provider = self.get_dataset_provider(path)
            dataset_info = provider.get_dataset_info(get_fs(path), path)
            self.cached_dataset_info = dataset_info
            self.cached_dataset_id = dataset_id
        return dataset_info

    def get_schema(self, dataset):
        path = dataset['url']
        provider = self.get_dataset_provider(path)
        schema_dict = provider.get_schema(get_fs(path), path)
        if 'summary' in dataset:
            schema_dict['summary'] = dataset['summary']
        if 'markers' in schema_dict:
            schema_dict['markers_read_only'] = schema_dict.pop('markers')
        return schema_dict

    def read_dataset(self, dataset, keys=[]):
        path = dataset['url']
        provider = self.get_dataset_provider(path)
        return provider.read_dataset(get_fs(path), path, keys=keys, dataset=dataset)

    def get_result(self, dataset, result_id):
        path = dataset['url']
        provider = self.get_dataset_provider(path)
        return provider.get_result(get_fs(path), path, dataset=dataset, result_id=result_id)
