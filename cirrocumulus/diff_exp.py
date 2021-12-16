import itertools

import numpy as np
import pandas as pd
import scipy.stats


def _power(X, power):
    return X ** power if isinstance(X, np.ndarray) else X.power(power)


def asarray(x):
    return x.toarray() if scipy.sparse.issparse(x) else x


class DE:

    def __init__(self, series: pd.Series, nfeatures: int, batch_size: int, get_batch_fn, base: float = None,
                 one_vs_rest: bool = True):
        """

        :param series: Categorical series in adata.obs to group by
        :param nfeatures: Number of features in adata
        :param batch_size: Number of features per batch
        :param get_batch_fn: Function to retrieve data from a batch
        :param base: adata.uns['log1p']['base']
        :param one_vs_rest: Whether to compare each group vs rest or all pairs of groups
        """
        mean_df = None
        variance_df = None
        frac_expressed_df = None
        indicator_df = pd.get_dummies(series)
        if one_vs_rest:
            pairs = []
            rest_indicator_df = pd.DataFrame()
            for c in indicator_df:
                rest_name = str(c) + '_rest'
                if rest_name in indicator_df:
                    counter = 1
                    rest_name = str(c) + '_rest-{}'.format(counter)
                    while rest_name in indicator_df:
                        counter = counter + 1
                        rest_name = str(c) + '_rest-{}'.format(counter)
                pairs.append((c, rest_name))
                rest_indicator_series = indicator_df[c].astype(bool)
                rest_indicator_series = ~rest_indicator_series
                rest_indicator_df[rest_name] = rest_indicator_series.astype(int)
            indicator_df = indicator_df.join(rest_indicator_df)
        else:
            pairs = list(itertools.combinations(series.cat.categories, 2))
        count_ = indicator_df.sum(axis=0)  # count per group
        A = scipy.sparse.coo_matrix(indicator_df.astype(float).T)
        n_row = A.shape[0]
        row_sums = np.asarray(A.sum(axis=1))
        D = scipy.sparse.dia_matrix(((row_sums.T ** -1), [0]), shape=(n_row, n_row))
        A = D * A
        dof = 1

        for i in range(0, nfeatures, batch_size):
            adata_batch = get_batch_fn(i)
            X = adata_batch.X
            mean_ = asarray(A @ X)  # (groups, genes)
            mean_sq = asarray(A @ _power(X, 2))
            sq_mean = mean_ ** 2
            var_ = mean_sq - sq_mean
            # enforce R convention (unbiased estimator) for variance
            precision = 2 << (42 if X.dtype == np.float64 else 20)
            # detects loss of precision in mean_sq - sq_mean, which suggests variance is 0
            var_[precision * var_ < sq_mean] = 0
            if dof != 0:
                var_ *= (count_ / (count_ - dof))[:, np.newaxis]

            frac_expressed_ = None
            if scipy.sparse.issparse(X):
                frac_expressed_ = asarray(A @ (X != 0))
            _mean_df = pd.DataFrame(mean_, columns=adata_batch.var.index, index=indicator_df.columns)
            _variance_df = pd.DataFrame(var_, columns=adata_batch.var.index, index=indicator_df.columns)

            # groups on rows, genes on columns
            mean_df = pd.concat((mean_df, _mean_df), axis=1) if mean_df is not None else _mean_df
            variance_df = pd.concat((variance_df, _variance_df), axis=1) if variance_df is not None else _variance_df
            if frac_expressed_ is not None:
                _frac_expressed_df = pd.DataFrame(frac_expressed_, columns=adata_batch.var.index,
                                                  index=indicator_df.columns)
                frac_expressed_df = pd.concat((frac_expressed_df, _frac_expressed_df),
                                              axis=1) if frac_expressed_df is not None else _frac_expressed_df

        if base is not None:
            expm1_func = lambda x: np.expm1(x * np.log(base))
        else:
            expm1_func = np.expm1
        pair2results = dict()
        for p in pairs:
            group_one, group_two = p
            nobs1 = count_.loc[group_one]
            nobs2 = count_.loc[group_two]

            # add small value to remove 0's
            foldchanges = np.log2(
                (expm1_func(mean_df.loc[group_one].values) + 1e-9) / (expm1_func(mean_df.loc[group_two].values) + 1e-9))
            with np.errstate(invalid="ignore"):
                scores, pvals = scipy.stats.ttest_ind_from_stats(
                    mean1=mean_df.loc[group_one],
                    std1=np.sqrt(variance_df.loc[group_one]),
                    nobs1=nobs1,
                    mean2=mean_df.loc[group_two],
                    std2=np.sqrt(variance_df.loc[group_two]),
                    nobs2=nobs2,
                    equal_var=False,  # Welch's
                )

            scores[np.isnan(scores)] = 0
            pvals[np.isnan(pvals)] = 1
            key = p[0] if one_vs_rest else p
            pair2results[key] = dict(scores=scores, pvals=pvals, logfoldchanges=foldchanges,
                                     frac_expressed1=frac_expressed_df.loc[
                                         group_one].values if frac_expressed_df is not None else None,
                                     frac_expressed2=frac_expressed_df.loc[
                                         group_two].values if frac_expressed_df is not None else None)
        self.pair2results = pair2results
