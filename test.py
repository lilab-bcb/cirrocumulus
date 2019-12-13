import anndata
import numpy as np
import scipy.sparse


X = scipy.sparse.csr_matrix(np.random.random((10, 4)))
data = anndata.AnnData(X=X)
data.write('test2.h5ad')
print(data.X.mean(axis=0).A1)

backed_data = anndata.read('test2.h5ad', backed='r')
print(backed_data.X[:, [0, 1]].mean(axis=0).A1)
