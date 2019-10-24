import pandas as pd

from .empty_cache import EmptyCache
from .file_system import FileSystem


class DatasetAPI:
    def __init__(self):
        self.suffix_to_provider = {}
        self.fs = FileSystem()
        self.caching_api = EmptyCache()

    def add(self, suffixes, provider):
        for suffix in suffixes:
            self.suffix_to_provider[suffix.lower()] = provider

    def schema(self, dataset):

        cached_schema = self.caching_api.get_schema(dataset.id)
        if cached_schema is not None:
            return cached_schema
        path = dataset['url']
        provider = self.suffix_to_provider[path[path.rfind('.') + 1:].lower()]
        value = provider.schema(self.fs, path)
        self.caching_api.cache_schema(dataset.id, value)
        return value

    def get_df(self, dataset, keys, embedding, index=False):
        path = dataset['url']
        provider = self.suffix_to_provider[path[path.rfind('.') + 1:].lower()]
        cached_dict = self.caching_api.get_df(dataset.id, keys, embedding, index=index)
        if keys is not None:
            keys = keys.copy()
        for column in cached_dict:
            if column in keys:
                keys.remove(column)
            elif column == 'index':
                index = False
        if embedding is not None:
            cached_embedding = True
            for column in embedding['coordinate_columns']:
                # check to see if all embedding columns are present
                if column not in cached_dict:
                    cached_embedding = False
                    break
            if cached_embedding:
                embedding = None

        df = provider.get_df(self.fs, path, keys, embedding, index=index)
        for column in cached_dict:
            if column not in df:
                df[column] = cached_dict[column]
        embedding_names = []
        if embedding is not None:
            for i in range(embedding['dimensions']):
                embedding_names.append(embedding['name'] + '_' + str(i + 1))
        for column in df:
            if not pd.api.types.is_numeric_dtype(df[column]) and not pd.api.types.is_categorical_dtype(df[column]):
                df[column] = df[column].astype('category')
        self.caching_api.cache_df(dataset.id, df.drop(cached_dict.keys(), axis=1))
        return df
