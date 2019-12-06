import pandas as pd


class IdsAggregator:

    def __init__(self):
        self.series = None

    def collect(self):
        return self.series

    def add(self, adata):
        # df index can be cell index or bin number
        # keep dataframe in long form
        series = adata.obs['id']
        self.series = pd.concat((self.series, series)) if self.series is not None else series
