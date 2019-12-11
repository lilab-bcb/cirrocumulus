class IdsAggregator:

    def execute(self, adata):
        return adata.obs['id'].to_list()
