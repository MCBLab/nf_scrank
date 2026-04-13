#!/usr/bin/env Rscript
library(Seurat)
library(dplyr)
library(scRank)

args <- commandArgs(trailingOnly = TRUE)

seuratObj <- args[1]
n_cores <- args[2]
species <- args[3]
target <- args[4]
column <- args[5]

target <- read.table(target, header = FALSE)
target <- target$V1

n_cores <- as.integer(n_cores)

cell_type <- sub(".RDS", "", seuratObj)

sc_obj <- readRDS(seuratObj)

obj <- CreateScRank(input = seuratObj,
                    species = species,
                    cell_type = column,
                    target = target[1])

obj <- Constr_net(obj)

weight <- obj@net[cell_type]
weight <- weight[colnames(weight),]

n_cells <- dim(sc_obj)[2]

# Save the object
saveRDS(weight, file = paste0(cell_type, "_weight_SCRANK_", n_cells, ".rds"))
