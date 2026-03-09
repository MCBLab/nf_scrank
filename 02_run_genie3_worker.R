#!/usr/bin/env Rscript
library(Seurat)
library(GENIE3)

args <- commandArgs(trailingOnly = TRUE)
rds_file <- args[1]
n_cores <- as.numeric(args[2])

message(">>> Processando arquivo: ", rds_file)
seuratObj <- readRDS(rds_file)
genes_use <- VariableFeatures(seuratObj)
target_gene <- "Smo"

if (target_gene %in% rownames(seuratObj) && !(target_gene %in% genes_use)) {
  genes_use <- c(genes_use, target_gene)
}

expr_mat <- as.matrix(GetAssayData(seuratObj, layer = "data")[genes_use, ])

message(">>> Iniciando GENIE3...")
weight_mat <- GENIE3(expr_mat, nCores = n_cores, verbose = TRUE)

out_name <- sub("_subset.rds$", "_weights.rds", rds_file)
if (out_name == rds_file) out_name <- sub(".rds$", "_weights.rds", rds_file)

saveRDS(weight_mat, file = out_name)
