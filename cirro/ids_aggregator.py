import pandas as pd


class IdsAggregator:

    def __init__(self):
        self.series = None

    def collect(self):
        return self.series

    def add(self, df):
        # df index can be cell index or bin number
        # keep dataframe in long form
        series = df[df['__selected']]['index']
        self.series = pd.concat((self.series, series)) if self.series is not None else series
