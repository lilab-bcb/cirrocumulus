from natsort import natsorted

from cirro.simple_data import SimpleData


def fraction_expressed(g):
    return (g > 0).sum() / len(g)


class DotPlotAggregator:

    def __init__(self, obs_measures, var_measures, dimensions):
        self.category_to_df = {}
        self.obs_measures = obs_measures
        self.var_measures = var_measures
        self.dimensions = dimensions

    def execute(self, adata):
        results = []
        # {categories:[], name:'', values:[{name:'', fractionExpressed:0, mean:0}]}
        adata_df = SimpleData.to_df(adata, self.obs_measures, self.var_measures, self.dimensions)

        for dimension in self.dimensions:
            df = adata_df.groupby(dimension).agg(['mean', fraction_expressed])
            df.index = df.index.astype('object')
            sorted_categories = natsorted(df.index)
            df = df.loc[sorted_categories]
            values = []
            dotplot_result = {'categories': df.index.values.tolist(), 'name': dimension, 'values': values}

            for measure in self.obs_measures + self.var_measures:
                series = df[measure]
                values.append({'name': measure, 'fractionExpressed': series['fraction_expressed'].to_list(),
                               'mean': series['mean'].to_list()})
            results.append(dotplot_result)
        return results
