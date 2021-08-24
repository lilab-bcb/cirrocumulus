import numpy as np
import pandas as pd
import scipy.stats

from cirrocumulus.groupby import GroupBy


class DE:

    def __init__(self, adata, obs_field, nfeatures, batch_size, get_batch_fn, key_set=None):
        group_by = GroupBy(adata, obs_field, key_set=key_set)
        A, keys = group_by.sparse_aggregator(True)
        mean = None
        variance = None
        count = None
        frac_expressed = None
        for i in range(0, nfeatures, batch_size):
            group_by.adata = get_batch_fn(i)  # hack to update anndata in groupby

            result = group_by.count_mean_var(A, keys)
            # groups on rows, genes on columns
            mean = pd.concat((mean, result['mean']), axis=1) if mean is not None else result['mean']
            variance = pd.concat((variance, result['var']), axis=1) if variance is not None else result['var']
            if count is None:
                count = result['count']
            if result['frac_expressed'] is not None:
                frac_expressed = pd.concat((frac_expressed, result['frac_expressed']),
                                           axis=1) if frac_expressed is not None else result['frac_expressed']
        mean1 = mean.iloc[0]
        mean2 = mean.iloc[1]
        if 'log1p' in adata.uns_keys() and adata.uns['log1p']['base'] is not None:
            expm1_func = lambda x: np.expm1(x * np.log(adata.uns['log1p']['base']))
        else:
            expm1_func = np.expm1
        # foldchanges = (expm1_func(mean1) + 1e-9) / (expm1_func(mean2) + 1e-9)  # add small value to remove 0's
        foldchanges = np.log2((expm1_func(mean1) + 1e-9) / (expm1_func(mean2) + 1e-9))
        variance1 = variance.iloc[0]
        variance2 = variance.iloc[1]
        nobs1 = count.iloc[0]
        nobs2 = count.iloc[1]
        with np.errstate(invalid="ignore"):
            scores, pvals = scipy.stats.ttest_ind_from_stats(
                mean1=mean1,
                std1=np.sqrt(variance1),
                nobs1=nobs1,
                mean2=mean2,
                std2=np.sqrt(variance2),
                nobs2=nobs2,
                equal_var=False,  # Welch's
            )

        scores[np.isnan(scores)] = 0
        pvals[np.isnan(pvals)] = 1
        self.scores = scores
        self.pvals = pvals
        self.logfoldchanges = foldchanges
        self.frac_expressed = frac_expressed
