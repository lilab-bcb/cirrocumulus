import pandas as pd

from .file_system import FileSystem


class DatasetAPI:
    def __init__(self):
        self.suffix_to_provider = {}
        self.fs = FileSystem()

    def add(self, suffixes, provider):
        for suffix in suffixes:
            self.suffix_to_provider[suffix.lower()] = provider

    def schema(self, path):
        provider = self.suffix_to_provider[path[path.rfind('.') + 1:].lower()]
        return provider.schema(self.fs, path)

    def get_df(self, path, keys, embedding):
        provider = self.suffix_to_provider[path[path.rfind('.') + 1:].lower()]
        df = provider.get_df(self.fs, path, keys, embedding)
        for column in df:
            if not pd.api.types.is_numeric_dtype(df[column]) and not pd.api.types.is_categorical_dtype(df[column]):
                df[column] = df[column].astype('category')
        return df
