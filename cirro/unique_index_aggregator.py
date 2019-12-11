class UniqueIndexAggregator:

    def execute(self, adata):
        # df index can be cell index or bin number
        # keep track of unique bins or indices
        count = adata.shape[0]
        index = adata.obs.index.unique()
        index.values.sort()
        return {'indices_or_bins': index.to_list(), 'count': count}
