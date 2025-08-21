import logging

import pandas as pd
import anndata
import scipy.sparse
from anndata import AnnData
from pandas import CategoricalDtype

from cirrocumulus.abstract_dataset import AbstractDataset
from cirrocumulus.anndata_util import ADATA_LAYERS_UNS_KEY, ADATA_MODULE_UNS_KEY, dataset_schema
from cirrocumulus.io_util import add_spatial, read_star_fusion_file


logger = logging.getLogger("cirro")

CATEGORICAL_FIELDS_CONVERT = ["seurat_clusters"]


def read_adata(path, filesystem, backed=False, spatial_directory=None, use_raw=False):
    path = path.rstrip(filesystem.sep)
    path_lc = path.lower()
    if path_lc.endswith(".loom"):
        adata = anndata.read_loom(filesystem.open(path))
    elif path_lc.endswith(".zarr"):
        adata = anndata.read_zarr(filesystem.get_mapper(path))
    elif path_lc.endswith(".tsv"):
        adata = read_star_fusion_file(filesystem.open(path))
    elif path_lc.endswith(".h5"):
        import scanpy as sc

        adata = sc.read_10x_h5(path)  # fsspec not supported
        sep_index = path.rfind(filesystem.sep)
        cells_path = "cells.parquet"
        analysis_path = "analysis"
        if sep_index != -1:
            cells_path = path[0:sep_index] + filesystem.sep + cells_path
            analysis_path = path[0:sep_index] + filesystem.sep + analysis_path

        if filesystem.exists(cells_path):  # xenium
            df = pd.read_parquet(cells_path).set_index("cell_id")
            df.index = df.index.astype(str)
            adata.obs = adata.obs.join(df)
            if "x_centroid" in adata.obs.columns:
                adata.obsm["X_spatial"] = adata.obs[["x_centroid", "y_centroid"]].values

                del adata.obs["x_centroid"]
                del adata.obs["y_centroid"]
        if filesystem.exists(analysis_path) and filesystem.isdir(analysis_path):  # xenium
            clustering_results = filesystem.glob(
                analysis_path + filesystem.sep + "clustering" + filesystem.sep + "**clusters.csv"
            )
            for clustering_result in clustering_results:
                df = pd.read_csv(clustering_result, index_col="Barcode")
                if len(df.columns) == 1:
                    df.index = df.index.astype(str)
                    df[df.columns[0]] = df[df.columns[0]].astype("category")
                    rename = dict()
                    # /analysis/clustering/gene_expression_graphclust/clusters.csv'
                    rename[df.columns[0]] = clustering_result.split(filesystem.sep)[-2]
                    df = df.rename(rename, axis=1)
                    adata.obs = adata.obs.join(df)
            umap_results = filesystem.glob(
                analysis_path + filesystem.sep + "umap" + filesystem.sep + "**projection.csv"
            )
            for umap_result in umap_results:
                df = pd.read_csv(umap_result, index_col="Barcode")
                if len(df.columns) == 2 or len(df.columns) == 3:
                    df.index = df.index.astype(str)
                    adata.obs = adata.obs.join(df)
                    adata.obsm["X_" + "_".join(umap_result.split(filesystem.sep)[-2:])] = adata.obs[
                        df.columns
                    ].values
                    for c in df.columns:
                        del adata.obs[c]

    elif path_lc.endswith(".rds"):  # Seurat, convert to h5ad
        h5_file = path + ".h5ad"
        import os

        if (
            not os.path.exists(h5_file)
            or abs(os.path.getmtime(h5_file) - os.path.getmtime(path)) > 0.00001
        ):
            import shutil
            import subprocess

            import pkg_resources

            print("Converting Seurat object")
            if os.path.exists(h5_file):
                os.remove(h5_file)
            subprocess.check_call(
                [
                    "Rscript",
                    pkg_resources.resource_filename("cirrocumulus", "seurat2h5ad.R"),
                    path,
                    h5_file,
                ]
            )
            shutil.copystat(path, h5_file)
        adata = anndata.read_h5ad(h5_file, backed="r" if backed else None)
        if adata.raw is not None and adata.shape[0] == adata.raw.shape[0]:
            print("Using adata.raw")
            adata = anndata.AnnData(
                X=adata.raw.X, var=adata.raw.var, obs=adata.obs, obsm=adata.obsm, uns=adata.uns
            )
    else:
        if backed:
            adata = anndata.read_h5ad(path, backed="r")
        else:
            adata = anndata.read_h5ad(filesystem.open(path))

    if "module" in adata.uns:
        adata.uns[ADATA_MODULE_UNS_KEY] = anndata.AnnData(
            X=adata.uns["module"]["X"], var=adata.uns["module"]["var"]
        )

    if "images" in adata.uns:
        images = adata.uns["images"]
        if isinstance(images, dict):
            adata.uns["images"] = [images]

    if use_raw and adata.raw is not None and adata.shape[0] == adata.raw.shape[0]:
        logger.info("Using adata.raw")
        adata = anndata.AnnData(
            X=adata.raw.X, var=adata.raw.var, obs=adata.obs, obsm=adata.obsm, uns=adata.uns
        )
    adata.var_names_make_unique()
    if spatial_directory is not None:
        if not add_spatial(adata, spatial_directory):
            logger.info("No spatial data found in {}".format(spatial_directory))

    for field in CATEGORICAL_FIELDS_CONVERT:
        if field in adata.obs and not isinstance(adata.obs[field].dtype, CategoricalDtype):
            logger.info("Converting {} to categorical".format(field))
            adata.obs[field] = adata.obs[field].astype(str).astype("category")
    return adata


class AnndataDataset(AbstractDataset):
    def __init__(self, backed=None):
        super().__init__()
        self.path_to_data = {}
        self.backed = backed

    def get_suffixes(self):
        return ["h5ad", "loom", "rds", "zarr", "h5"]

    def get_result(self, filesystem, path, dataset, result_id):
        adata = self.get_data(filesystem, path)
        if result_id in adata.uns:
            return str(adata.uns[result_id])
        return super().get_result(filesystem, path, dataset, result_id)

    def read_adata(self, filesystem, path):
        adata = read_adata(path, filesystem, self.backed)
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
        if "cirro-schema" in adata.uns:
            import json

            s = json.loads(str(adata.uns["cirro-schema"]))
            keys = ["markers", "results"]
            for key in keys:
                if key in s:
                    schema[key] = s[key]
        return schema

    @staticmethod
    def get_X(adata, keys, layer=None):
        if len(keys) == 1 and isinstance(keys[0], slice):  # special case if slice specified
            keys = keys[0]
        d = adata[:, keys]
        X = d.X if layer is None else d.layers[layer]
        if scipy.sparse.issparse(X) and not scipy.sparse.isspmatrix_csc(X):
            X = X.tocsc()

        var = pd.DataFrame(index=d.var.index)
        return X, var

    def read_dataset(self, filesystem, path, keys=None, dataset=None):
        adata = self.get_data(filesystem, path)
        if keys is None:
            keys = {}
        keys = keys.copy()
        X_keys = keys.pop("X", [])
        obs_keys = keys.pop("obs", [])
        basis_keys = keys.pop("basis", [])
        module_keys = keys.pop("module", [])
        X = None
        obs = None
        var = None
        obsm = {}
        adata_modules = None
        layers = {}
        for layer_key in keys.keys():
            X_layer, var_layer = AnndataDataset.get_X(adata, keys[layer_key], layer_key)
            adata_layer = AnnData(X=X_layer, var=var_layer)
            layers[layer_key] = adata_layer
        if len(X_keys) > 0:
            X, var = AnndataDataset.get_X(adata, X_keys)
        if len(module_keys) > 0:
            if len(module_keys) == 1 and isinstance(
                module_keys[0], slice
            ):  # special case if slice specified
                module_keys = module_keys[0]
            adata_modules = adata.uns[ADATA_MODULE_UNS_KEY][:, module_keys]

        if len(obs_keys) > 0:
            obs = pd.DataFrame()
            for key in obs_keys:
                if key == "index":
                    values = adata.obs.index.values
                else:
                    values = adata.obs[key].values
                obs[key] = values
        if len(basis_keys) > 0:
            for key in basis_keys:
                embedding_data = adata.obsm[key]
                obsm[key] = embedding_data

        if X is None and obs is None and len(obsm.keys()) == 0:
            obs = pd.DataFrame(index=pd.RangeIndex(adata.shape[0]).astype(str))
        adata = AnnData(X=X, obs=obs, var=var, obsm=obsm)
        if adata_modules is not None:
            adata.uns[ADATA_MODULE_UNS_KEY] = adata_modules
        adata.uns[ADATA_LAYERS_UNS_KEY] = layers
        return adata
