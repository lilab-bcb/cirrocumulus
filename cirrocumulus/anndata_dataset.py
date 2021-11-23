import anndata
import pandas as pd
import scipy.sparse
from anndata import AnnData

from cirrocumulus.abstract_dataset import AbstractDataset
from cirrocumulus.anndata_util import dataset_schema, ADATA_MODULE_UNS_KEY
from cirrocumulus.io_util import read_star_fusion_file


class AnndataDataset(AbstractDataset):

    def __init__(self, backed=None):
        super().__init__()
        self.path_to_data = {}
        self.backed = backed

    def get_suffixes(self):
        return ['h5ad', 'loom', 'rds', 'zarr']

    def get_result(self, filesystem, path, dataset, result_id):
        adata = self.get_data(filesystem, path)
        if result_id in adata.uns:
            return str(adata.uns[result_id])
        return super().get_result(filesystem, path, dataset, result_id)

    def read_adata(self, filesystem, path):
        path_lc = path.lower()
        path_lc = path_lc.rstrip('/')
        if path_lc.endswith('.loom'):
            adata = anndata.read_loom(filesystem.open(path))
        elif path_lc.endswith('.zarr'):
            adata = anndata.read_zarr(filesystem.get_mapper(path))
        elif path_lc.endswith('.tsv'):
            adata = read_star_fusion_file(filesystem.open(path))
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
        else:
            adata = anndata.read_h5ad(filesystem.open(path), backed='r' if self.backed else None)
        if 'module' in adata.uns:
            adata.uns[ADATA_MODULE_UNS_KEY] = anndata.AnnData(X=adata.uns['module']['X'],
                                                              var=adata.uns['module']['var'])
        return adata

    def add_data(self, path, data):
        self.path_to_data[path] = data

    def get_data(self, filesystem, path):
        adata = self.path_to_data.get(path)
        if adata is None:
            adata = self.read_adata(filesystem, path)
            if scipy.sparse.isspmatrix_csr(adata.X) and adata.X.shape[1] > 1:
                adata.X = adata.X.tocsc()
            self.add_data(path, adata)
        return adata

    def get_schema(self, filesystem, path):
        adata = self.get_data(filesystem, path)
        schema = dataset_schema(adata)
        if 'cirro-schema' in adata.uns:
            import json
            s = json.loads(str(adata.uns['cirro-schema']))
            keys = ['markers', 'results']
            for key in keys:
                if key in s:
                    schema[key] = s[key]
        return schema

    def read_dataset(self, filesystem, path, keys=None, dataset=None):
        adata = self.get_data(filesystem, path)
        if keys is None:
            keys = {}
        keys = keys.copy()
        X_keys = keys.pop('X', [])
        obs_keys = keys.pop('obs', [])
        basis_keys = keys.pop('basis', [])
        module_keys = keys.pop('module', [])
        X = None
        obs = None
        var = None
        obsm = {}
        adata_modules = None
        if len(X_keys) > 0:
            if len(X_keys) == 1 and isinstance(X_keys[0], slice):  # special case if slice specified
                X_keys = X_keys[0]
            d = adata[:, X_keys]
            if scipy.sparse.issparse(d.X) and not scipy.sparse.isspmatrix_csc(d.X):
                d.X = d.X.tocsc()
            X = d.X
            var = pd.DataFrame(index=d.var.index)
        if len(module_keys) > 0:
            if len(module_keys) == 1 and isinstance(module_keys[0], slice):  # special case if slice specified
                module_keys = module_keys[0]
            adata_modules = adata.uns[ADATA_MODULE_UNS_KEY][:, module_keys]

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
            if X is None:  # anndata requires empty X
                X = scipy.sparse.coo_matrix(([], ([], [])), shape=(embedding_data.shape[0], 0))
        if X is None and obs is None and len(obsm.keys()) == 0:
            obs = pd.DataFrame(index=pd.RangeIndex(adata.shape[0]).astype(str))

        adata = AnnData(X=X, obs=obs, var=var, obsm=obsm)
        if adata_modules is not None:
            adata.uns[ADATA_MODULE_UNS_KEY] = adata_modules
        return adata
