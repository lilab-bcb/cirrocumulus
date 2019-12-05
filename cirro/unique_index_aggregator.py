class UniqueIndexAggregator:

    def __init__(self):
        self.index = None
        self.count = 0

    def collect(self):
        self.index.values.sort()
        return {'indices_or_bins': self.index.values.tolist(), 'count': self.count}

    def add(self, adata):
        # df index can be cell index or bin number
        # keep track of unique bins or indices
        self.count += adata.shape[0]
        self.index = self.index.union(adata.obs.index.unique()) if self.index is not None else adata.obs.index.unique()
