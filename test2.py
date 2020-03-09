import dask.dataframe as dd
import pandas as pd
import scipy.sparse

c = ['a', 'b', 'a', 'c', 'b', 'a', 'a', 'b', 'e', 'x']
test_df = pd.DataFrame(scipy.sparse.eye(10).toarray())
test_df['x'] = c
test_df['x'] = test_df['x'].astype('category')
test_df['y'] = 1
test_dd = dd.from_pandas(test_df, 2)

foo = None


def mode_chunk(x):
    r = x.value_counts(sort=False)
    global foo
    foo = r
    return r


def mode_agg(x):
    print("XXXXX")
    return x.apply(lambda s: s.groupby(level=-1).sum())


def mode_finalize(x):
    # return x.groupby(level=level).apply(lambda s: s.reset_index(level=level, drop=True).idxmax())
    return x.groupby(level=0).agg(lambda s: s.nlargest(1).index[0][1])


dask_agg_mode = dd.Aggregation(
    name='dd_mode',
    chunk=mode_chunk,
    agg=mode_agg,
    finalize=mode_finalize
)
test_dd.groupby(0).agg({'x': dask_agg_mode, 1: 'min'}).compute()
test_dd.groupby(0).agg({1: 'min'}).compute()  # works
g = test_dd.groupby(0)
g.agg({'x': dask_agg_mode}).compute()  # fails

grouped = df.groupby('batch')
X = adata.X
X = X[:, [0, 1]]
for key, g in grouped:
    indices = grouped.indices[key]
    X_group = X[indices]
    X_summary = X_group.mean(axis=0).A1
