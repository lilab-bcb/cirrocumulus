class IdsAggregator:

    def execute(self, adata):
        return adata.obs['index']
