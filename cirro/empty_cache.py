class EmptyCache:

    def cache_schema(self, schema, dataset_id):
        pass

    def get_schema(self, dataset_id):
        return None

    def cache_df(self, df, dataset_id):
        pass

    def get_df(self, dataset_id, keys, embedding_key=None, index=False):
        return {}
