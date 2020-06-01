import pandas as pd


class UniqueAggregator:

    def __init__(self, column):
        self.column = column

    def execute(self, adata):
        # column can be cell index or bin number
        if self.column == 'index':
            values = adata.obs.index
        else:
            values = adata.obs[self.column].unique()
            if isinstance(values, pd.arrays.SparseArray):
                values = values.to_dense()
            values.sort()
        return {'indices_or_bins': values}
