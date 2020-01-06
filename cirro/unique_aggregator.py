class UniqueAggregator:

    def __init__(self, column):
        self.column = column

    def execute(self, adata):
        # column can be cell index or bin number
        if self.column == 'index':
            values = adata.obs.index.to_list()
        else:
            values = adata.obs[self.column].unique()
            values.sort()
            values = values.tolist()
        return {'indices_or_bins': values}
