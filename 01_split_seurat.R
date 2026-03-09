#!/usr/bin/env Rscript
library(Seurat)
args <- commandArgs(trailingOnly = TRUE)
seurat_file <- args[1]

message(">>> Lendo: ", seurat_file)
seuratObj <- readRDS(seurat_file)

message(">>> Filtrando Clusters 0 e 1...")
Idents(seuratObj) <- "seurat_clusters"
seuratObj <- subset(seuratObj, idents = c("0", "1"))
seuratObj <- RenameIdents(seuratObj, c("1" = "NodeB_Sensivel", "0" = "NodeA_Resistente"))

# 2000 genes variáveis para o GENIE3 ter bastante informação
seuratObj <- FindVariableFeatures(seuratObj, nfeatures = 2000)

dir.create("nextflow_inputs", showWarnings = FALSE)

for (ct in unique(Idents(seuratObj))) {
  message("   -> Exportando: ", ct)
  saveRDS(subset(seuratObj, idents = ct), file = paste0("nextflow_inputs/", ct, "_subset.rds"))
}
message(">>> Divisão concluída com sucesso!")
