#!/usr/bin/env Rscript

method = Sys.getenv('CIRRO_SEURAT_CONVERTER') # either reticulate or SeuratDisk
if (method == '') {
  method <- 'reticulate'
}

# if (method == 'reticulate' && !requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# if (method == 'reticulate' && !require('SingleCellExperiment', quietly = TRUE))
#   BiocManager::install("SingleCellExperiment")
if (!require('Seurat', quietly = TRUE))
  install.packages('Seurat', repos = 'https://cloud.r-project.org/')

if (method == 'SeuratDisk' && !require('remotes', quietly = TRUE))
  install.packages('remotes', repos = 'https://cloud.r-project.org/')

if (method == 'SeuratDisk' && !require('SeuratDisk', quietly = TRUE))
  remotes::install_github("mojaveazure/seurat-disk")

library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
rds <- readRDS(args[1])
rds <- UpdateSeuratObject(rds)
h5ad_path <- args[2]
if (file.exists(h5ad_path)) {
  unlink(h5ad_path)
}

if (method == 'reticulate') {
  library(reticulate)
  library(Matrix)
  library(ps)

  anndata <- import("anndata")
  assay <- DefaultAssay(object = rds)
  exprs <- GetAssayData(object = rds, slot = "data", assay = assay)
  col_data <- rds[[]]
  col_data$ident <- Idents(object = rds)
  row_data <- rds[[assay]][[]]
  if (length(row_data) == 0) {
    row_data['tmp'] = 'tmp'
  }
  adata <- anndata$AnnData(X = t(exprs), obs = col_data, var = row_data)
  # dim_names = reducedDimNames(sce)
  # for (dim_name in dim_names) {
  #   adata$obsm$setdefault(dim_name, reducedDim(sce, dim_name))
  # }
  for (dr in Seurat:::FilterObjects(object = rds, classes.keep = "DimReduc")) {
    adata$obsm$setdefault(dr, Embeddings(object = rds[[dr]]))
  }
  adata$write(h5ad_path)

  # Fix hanging R process: There appear to be 1 leaked semaphore objects to clean up at shutdown.
  child_processes <- ps_children()
  if (length(child_processes) == 1) {
    ps_kill(child_processes[[1]])
  }
} else {
  #if (!require('loomR', quietly = TRUE))
  #  remotes::install_github(repo = 'mojaveazure/loomR', ref = 'develop')
  #library(loomR)
  library(SeuratDisk)
  #loom <- as.loom(rds, filename = args[2], verbose = FALSE)
  #loom$close_all()
  h5_seurat_path <- paste0(tools::file_path_sans_ext(h5ad_path), '.h5Seurat')
  SaveH5Seurat(rds, filename = h5_seurat_path)
  Convert(h5_seurat_path, dest = "h5ad")
  unlink(h5_seurat_path)
}


