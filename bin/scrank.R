#!/usr/bin/env Rscript
library(Seurat)
library(dplyr)
library(scRank)

args <- commandArgs(trailingOnly = TRUE)

seuratObj <- args[1]
species <- args[2]
targets <- args[3]
column <- args[4]
n_cores <- args[5]

targets <- readLines(targets)
target <- targets[1]

n_cores <- as.integer(n_cores)

cell_type <- sub(".RDS", "", seuratObj)

sc_obj <- readRDS(seuratObj)

obj <- CreateScRank(input = sc_obj,
                    species = species,
                    cell_type = column,
                    target = target)

obj <- Constr_net(obj)

weight <- obj@net[cell_type]
weight <- weight[colnames(weight),]

n_cells <- dim(sc_obj)[2]

# Save the object
saveRDS(weight, file = paste0(cell_type, "_weight_SCRANK_", n_cells, ".rds"))
