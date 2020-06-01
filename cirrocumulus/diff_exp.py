import numpy as np


def _ecdf(x):
    '''no frills empirical cdf used in fdrcorrection
    '''
    nobs = len(x)
    return np.arange(1, nobs + 1) / float(nobs)


def fdrcorrection(pvals, alpha=0.05, method='indep', is_sorted=False):
    '''pvalue correction for false discovery rate

    This covers Benjamini/Hochberg for independent or positively correlated and
    Benjamini/Yekutieli for general or negatively correlated tests. Both are
    available in the function multipletests, as method=`fdr_bh`, resp. `fdr_by`.

    Parameters
    ----------
    pvals : array_like
        set of p-values of the individual tests.
    alpha : float
        error rate
    method : {'indep', 'negcorr')

    Returns
    -------
    rejected : ndarray, bool
        True if a hypothesis is rejected, False if not
    pvalue-corrected : ndarray
        pvalues adjusted for multiple hypothesis testing to limit FDR

    Notes
    -----

    If there is prior information on the fraction of true hypothesis, then alpha
    should be set to alpha * m/m_0 where m is the number of tests,
    given by the p-values, and m_0 is an estimate of the true hypothesis.
    (see Benjamini, Krieger and Yekuteli)

    The two-step method of Benjamini, Krieger and Yekutiel that estimates the number
    of false hypotheses will be available (soon).

    Method names can be abbreviated to first letter, 'i' or 'p' for fdr_bh and 'n' for
    fdr_by.



    '''
    pvals = np.asarray(pvals)

    if not is_sorted:
        pvals_sortind = np.argsort(pvals)
        pvals_sorted = np.take(pvals, pvals_sortind)
    else:
        pvals_sorted = pvals  # alias

    if method in ['i', 'indep', 'p', 'poscorr']:
        ecdffactor = _ecdf(pvals_sorted)
    elif method in ['n', 'negcorr']:
        cm = np.sum(1. / np.arange(1, len(pvals_sorted) + 1))  # corrected this
        ecdffactor = _ecdf(pvals_sorted) / cm
    ##    elif method in ['n', 'negcorr']:
    ##        cm = np.sum(np.arange(len(pvals)))
    ##        ecdffactor = ecdf(pvals_sorted)/cm
    else:
        raise ValueError('only indep and negcorr implemented')

    pvals_corrected_raw = pvals_sorted / ecdffactor
    pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
    del pvals_corrected_raw
    pvals_corrected[pvals_corrected > 1] = 1
    if not is_sorted:
        pvals_corrected_ = np.empty_like(pvals_corrected)
        pvals_corrected_[pvals_sortind] = pvals_corrected
        del pvals_corrected

        return pvals_corrected_
    else:
        return pvals_corrected


def diff_exp(X, mask):
    pvals = np.full(X.shape[1], 1.0)
    tscores = np.full(X.shape[1], 0)
    mat_cond1 = X[mask]
    mat_cond2 = X[~mask]
    n1 = mat_cond1.shape[0]
    n2 = mat_cond2.shape[0]
    mean1 = mat_cond1.mean(axis=0).A1
    mean2 = mat_cond2.mean(axis=0).A1
    psum1 = mat_cond1.power(2).sum(axis=0).A1
    s1sqr = (psum1 - n1 * (mean1 ** 2)) / (n1 - 1)

    psum2 = mat_cond2.power(2).sum(axis=0).A1
    s2sqr = (psum2 - n2 * (mean2 ** 2)) / (n2 - 1)

    percent1 = (mat_cond1.getnnz(axis=0) / n1 * 100.0).astype(np.float32)

    percent2 = (mat_cond2.getnnz(axis=0) / n2 * 100.0).astype(np.float32)

    import scipy.stats as ss

    var_est = s1sqr / n1 + s2sqr / n2
    idx = var_est > 0.0
    if idx.sum() > 0:
        tscore = (mean1[idx] - mean2[idx]) / np.sqrt(var_est[idx])
        v = (var_est[idx] ** 2) / (
                (s1sqr[idx] / n1) ** 2 / (n1 - 1) + (s2sqr[idx] / n2) ** 2 / (n2 - 1)
        )
        pvals[idx] = ss.t.sf(np.fabs(tscore), v) * 2.0  # two-sided
        tscores[idx] = tscore
    #  qvals = fdrcorrection(pvals)
    log_fold_change = mean1 - mean2
    x_avg = (mean1 + mean2) / 2
    x_max = x_avg.max()
    x_min = x_avg.min() - 0.001  # to avoid divide by zero
    weights = (x_avg - x_min) / (x_max - x_min)
    WAD = log_fold_change * weights

    return dict(WAD=WAD, mean1=mean1, mean2=mean2, percent1=percent1, percent2=percent2, tscore=tscore)
