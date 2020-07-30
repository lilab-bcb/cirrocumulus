import pandas as pd


class UniqueAggregator:

    def __init__(self, column):
        self.column = column

    def execute(self, df):
        # column can be cell index or bin number
        if self.column == 'index':
            values = df.index
        else:
            values = df[self.column].unique()
            if isinstance(values, pd.arrays.SparseArray):
                values = values.to_dense()
            values.sort()
        return {'indices_or_bins': values}
