#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(WGCNA)
})

args <- commandArgs(trailingOnly = TRUE)
input_rds <- args[1]
column <- args[2]
n_cores <- as.numeric(args[3])
targets_file <- args[4]

allowWGCNAThreads(nThreads = n_cores)

cell_type <- sub("\\.RDS$|\\.rds$", "", basename(input_rds))
seuratObj <- readRDS(input_rds)
n_cells <- ncol(seuratObj)

if (n_cells < 150) {
  quit(save = "no", status = 0)
}

genes_4_use <- seuratObj@misc$gene4use

if (file.exists(targets_file)) {
  targets <- readLines(targets_file)
  targets <- targets[targets != ""]
  genes_4_use <- unique(c(genes_4_use, intersect(targets, rownames(seuratObj))))
}

seuratObj$wgcna_group <- "target_clone"
nome_arquivo <- paste0(cell_type, "_weight_hdWGCNA_", n_cells, ".rds")

seuratObj <- SetupForWGCNA(
  seuratObj,
  gene_select = "custom",
  features = genes_4_use,
  wgcna_name = "network"
)

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
power_est <- GetActiveWGCNA(seuratObj)$sft$powerEstimate

final_power <- if(is.null(power_est) || is.na(power_est) || is.infinite(power_est)) 9 else power_est

seuratObj <- ConstructNetwork(
  seuratObj, 
  setDatExpr = FALSE, 
  overwrite_tom = TRUE, 
  tom_name = "network", 
  soft_power = final_power,
  minModuleSize = 10 
)

weight <- as.matrix(GetTOM(seuratObj))
weight <- weight[colnames(weight), ]

saveRDS(weight, file = nome_arquivo)