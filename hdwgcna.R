suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(WGCNA)
})

allowWGCNAThreads(nThreads = 4)
file_in <- "C:/Users/leoga/projeto_nextflow/data/vehicle_final_server.rds"
out_dir <- "C:/Users/leoga/projeto_nextflow/data/"

seuratObj <- readRDS(file_in)
Idents(seuratObj) <- "clone_annotation"

seuratObj <- SetupForWGCNA(
  seuratObj,
  gene_select = "fraction", 
  fraction = 0.05,          
  wgcna_name = "Medullo_SHH"
)

seuratObj <- MetacellsByGroups(
  seurat_obj = seuratObj,
  group.by = c("clone_annotation", "orig.ident"),
  k = 25, 
  max_shared = 10, 
  ident.group = "clone_annotation"
)

seuratObj <- NormalizeMetacells(seuratObj)

valid_celltypes <- as.character(unique(seuratObj@misc$Medullo_SHH$wgcna_metacell_obj@meta.data$clone_annotation))

for (clone in valid_celltypes) {
  
  seuratObj <- SetDatExpr(
    seuratObj,
    group_name = clone, 
    group.by = "clone_annotation",
    assay = 'RNA',
    slot = 'data'
  )
  
  seuratObj <- TestSoftPowers(
    seuratObj,
    networkType = 'signed' 
  )
  
  seuratObj <- ConstructNetwork(
    seuratObj,
    setDatExpr = FALSE, 
    tom_name = paste0('TOM_', clone), 
    overwrite_tom = TRUE 
  )
  
  weight <- as.matrix(GetTOM(seuratObj))
  
  weight <- weight[colnames(weight), ]
  
  n_cells <- sum(seuratObj$clone_annotation == clone)

  nome_arquivo <- paste0(out_dir, clone, "_weight_hdWGCNA_", n_cells, ".rds")
  saveRDS(weight, file = nome_arquivo)
}