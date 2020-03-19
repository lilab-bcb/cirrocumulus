def mode_and_purity(x):
    value_counts = x.value_counts(sort=False)
    largest = value_counts.nlargest(1)
    purity = largest[0] / value_counts.sum()
    return largest.index[0], purity


def sparse_density(x):
    return x.sparse.density


def np_density(x):
    return (x != 0).sum() / len(x)
