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

mat <- as.matrix(sc_obj[sc_obj@misc$gene4use]@assays$RNA$counts)
mat <- mat[,colSums(mat) > 30]

n_cells <- dim(sc_obj)[2]

outputH0 <- scTenifoldNet(X = mat, Y = mat, nc_nNet = 10, nc_nCells = n_cells/2, nCores = n_cores, qc = F)
weight <- as.matrix(outputH0$tensorNetworks$X)

# Save the object
saveRDS(weight, file = paste0(cell_type, "_weight_sctenifoldnet_", n_cells, ".rds"))
