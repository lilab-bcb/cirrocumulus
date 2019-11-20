class UniqueIndexAggregator:

    def __init__(self):
        self.index = None
        self.count = 0

    def collect(self):
        return {'indices_or_bins': self.index.values.tolist(), 'count': self.count}

    def add(self, df):
        # df index can be cell index or bin number
        # keep track of unique bins or indices
        self.count += len(df)
        self.index = self.index.union(df.index.unique()) if self.index is not None else df.index.unique()
