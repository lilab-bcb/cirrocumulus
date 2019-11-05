class SelectionAggregator:

    def __init__(self):
        self.index = None
        self.count = 0

    def collect(self):
        return self.index, self.count


    def add(self, df):
        # df index can be cell index or bin number
        # keep track of unique bins or indices
        df = df[df['__selected']]
        self.count += len(df)
        self.index = self.index.union(df.index.unique()) if self.index is not None else df.index.unique()
