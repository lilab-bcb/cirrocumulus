#!/usr/bin/env Rscript

library(Seurat)
args <- commandArgs(trailingOnly = TRUE)
rds <- readRDS(args[1])
rds <- UpdateSeuratObject(rds)
h5ad_path <- args[2]
if (file.exists(h5ad_path)) {
  unlink(h5ad_path)
}

if (!require('remotes', quietly = TRUE))
  install.packages('remotes', repos = 'https://cloud.r-project.org/')

if (!require('SeuratDisk', quietly = TRUE))
  remotes::install_github("mojaveazure/seurat-disk")

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
