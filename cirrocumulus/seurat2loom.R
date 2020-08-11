#!/usr/bin/env Rscript

library(Seurat)

if (!require('remotes', quietly = TRUE))
  install.packages('remotes', repos = 'https://cloud.r-project.org/')
if (!require('loomR', quietly = TRUE))
  remotes::install_github(repo = 'mojaveazure/loomR', ref = 'develop')
library(loomR)
args <- commandArgs(trailingOnly = TRUE)
s <- readRDS(args[1])
if (file.exists(args[2])) {
  unlink(args[2])
}
loom <- as.loom(s, filename = args[2], verbose = FALSE)
loom$close_all()