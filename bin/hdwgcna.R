#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(WGCNA)
})

args <- commandArgs(trailingOnly = TRUE)
input_rds <- args[1]
n_cores <- as.numeric(args[2])
targets_file <- args[3]

allowWGCNAThreads(nThreads = n_cores)

cell_type <- sub("\\.RDS$|\\.rds$", "", basename(input_rds))
seuratObj <- readRDS(input_rds)
n_cells <- ncol(seuratObj)

if (n_cells < 150) {
  quit(save = "no", status = 0)
}

if (file.exists(targets_file)) {
  targets <- readLines(targets_file)
  targets <- targets[targets != ""]
  
  if (length(VariableFeatures(seuratObj)) == 0) {
    seuratObj <- FindVariableFeatures(seuratObj, selection.method = "vst", nfeatures = 2000)
  }
  
  current_hvgs <- VariableFeatures(seuratObj)
  valid_targets <- intersect(targets, rownames(seuratObj))
  VariableFeatures(seuratObj) <- unique(c(current_hvgs, valid_targets))
}

seuratObj$wgcna_group <- "target_clone"
nome_arquivo <- paste0(cell_type, "_weight_hdWGCNA_", n_cells, ".rds")

seuratObj <- SetupForWGCNA(seuratObj, wgcna_name = "network")
seuratObj <- MetacellsByGroups(
  seurat_obj = seuratObj,
  group.by = "wgcna_group",
  k = 25,
  min_cells = 100, 
  max_shared = 10,
  ident.group = "wgcna_group"
)
seuratObj <- NormalizeMetacells(seuratObj)

seuratObj <- SetDatExpr(
  seuratObj,
  group_name = "target_clone",
  group.by = "wgcna_group",
  assay = 'RNA',
  layer = 'data' 
)

seuratObj <- TestSoftPowers(seuratObj, networkType = 'signed')

net_data <- GetActiveWGCNA(seuratObj)
power_est <- net_data$sft$powerEstimate

if(is.null(power_est) || is.na(power_est) || is.infinite(power_est)) {
    seuratObj <- ConstructNetwork(seuratObj, setDatExpr = FALSE, overwrite_tom = TRUE, tom_name = "network", soft_power = 9)
} else {
    seuratObj <- ConstructNetwork(seuratObj, setDatExpr = FALSE, overwrite_tom = TRUE, tom_name = "network")
}

weight <- as.matrix(GetTOM(seuratObj))
weight <- weight[colnames(weight), ]
saveRDS(weight, file = nome_arquivo)
