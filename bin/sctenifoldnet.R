#!/usr/bin/env Rscript
library(Seurat)
library(dplyr)
library(scTenifoldNet)

args <- commandArgs(trailingOnly = TRUE)

seuratObj <- args[1]
n_cores <- args[2]
n_cores <- as.integer(n_cores)

cell_type <- sub(".RDS", "", seuratObj)

sc_obj <- readRDS(seuratObj)

mat <- as.matrix(sc_obj[sc_obj@misc$gene4use]@assays$RNA$data)

outputH0 <- scTenifoldNet(X = mat, Y = mat, nc_nNet = 4, nc_nCells = 10, nCores = n_cores, qc = F)
weight <- as.matrix(outputH0$tensorNetworks$X)

n_cells <- dim(sc_obj)[2]

# Save the object
saveRDS(weight, file = paste0(cell_type, "_weight_sctenifoldnet_", n_cells, ".rds"))
