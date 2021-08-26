import anndata
import pandas as pd
import scipy.sparse
from anndata import AnnData

from cirrocumulus.abstract_dataset import AbstractDataset
# only works with local files
from cirrocumulus.anndata_util import datasets_schema
from cirrocumulus.io_util import read_star_fusion_file


class AnndataDataset(AbstractDataset):

    def __init__(self, backed=None, force_sparse=True, extensions=['h5ad', 'loom', 'zarr', 'rds']):
        super().__init__()
        self.path_to_data = {}
        self.backed = backed
        self.force_sparse = force_sparse
        self.extensions = extensions

    def get_suffixes(self):
        return self.extensions

    def read_adata(self, path):
        path_lc = path.lower()
        if path_lc.endswith('.loom'):
            return anndata.read_loom(path)
        elif path_lc.endswith('.zarr'):
            return anndata.read_zarr(path)
        elif path_lc.endswith('.tsv'):
            return read_star_fusion_file(path)
        elif path_lc.endswith('.rds'):  # Seurat, convert to h5ad
            h5_file = path + '.h5ad'
            import os
            if not os.path.exists(h5_file) or abs(os.path.getmtime(h5_file) - os.path.getmtime(path)) > 0.00001:
                import subprocess
                import pkg_resources
                import shutil
                print('Converting Seurat object')
                if os.path.exists(h5_file):
                    os.remove(h5_file)
                subprocess.check_call(
                    ['Rscript', pkg_resources.resource_filename("cirrocumulus", 'seurat2h5ad.R'), path, h5_file])
                shutil.copystat(path, h5_file)
            adata = anndata.read_h5ad(h5_file, backed='r' if self.backed else None)
            if adata.raw is not None and adata.shape[0] == adata.raw.shape[0]:
                print('Using adata.raw')
                adata = anndata.AnnData(X=adata.raw.X, var=adata.raw.var, obs=adata.obs, obsm=adata.obsm, uns=adata.uns)
            return adata
        return anndata.read_h5ad(path, backed='r' if self.backed else None)
        # elif path.endswith('.mtx'):
        #
        #     return anndata.read_mtx(path, backed=self.backed)
        # elif path.endswith('.txt'):
        #
        #     return anndata.read_text(path, backed=self.backed)
        # elif path.endswith('.csv'):
        #
        #     return anndata.read_csv(path, backed=self.backed)

    def add_data(self, path, data):
        self.path_to_data[path] = data

    def get_data(self, path):
        adata = self.path_to_data.get(path)
        if adata is None:
            adata = self.read_adata(path)
            if scipy.sparse.isspmatrix_csr(adata.X) and adata.X.shape[1] > 1:
                adata.X = adata.X.tocsc()
            self.add_data(path, adata)
        return adata

    def schema(self, filesystem, path):
        return datasets_schema([self.get_data(path)])

    def read_dataset(self, filesystem, path, keys=None, dataset=None, schema=None):
        adata = self.get_data(path)
        if keys is None:
            keys = {}
        keys = keys.copy()
        var_keys = keys.pop('X', [])
        obs_keys = keys.pop('obs', [])
        basis_keys = keys.pop('basis', [])
        X = None
        obs = None
        var = None
        obsm = {}
        if len(var_keys) > 0:
            X = adata[:, var_keys].X
            if not scipy.sparse.issparse(X):
                X = scipy.sparse.csc_matrix(X)
            var = pd.DataFrame(index=var_keys)
        # for key in keys.keys():
        #     if df is None:
        #         df = pd.DataFrame()
        #     d = adata.uns[key]
        #     features = keys[key]
        #     X = d['X'][:, d['var'].index.get_indexer_for(features)]
        #     for i in range(len(features)):
        #         df[features[i]] = X[:, i]

        if len(obs_keys) > 0:
            obs = pd.DataFrame()
            for key in obs_keys:
                if key == 'index':
                    values = adata.obs.index.values
                else:
                    values = adata.obs[key].values
                obs[key] = values
        if len(basis_keys) > 0:
            for key in basis_keys:
                embedding_data = adata.obsm[key]
                obsm[key] = embedding_data
                if X is None:
                    X = scipy.sparse.coo_matrix(([], ([], [])), shape=(embedding_data.shape[0], 0))
        return AnnData(X=X, obs=obs, var=var, obsm=obsm)
