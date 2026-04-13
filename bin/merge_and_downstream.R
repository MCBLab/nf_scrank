#!/usr/bin/env Rscript
library(Seurat)
library(scRank)
library(dplyr)
library(GENIE3)

args <- commandArgs(trailingOnly = TRUE)

seuratObj_path <- args[1]
targets_file <- args[2]
species <- args[3]
column <- args[4]
rds_files <- args[5:length(args)]

cell_types <- sub("_weight.*", "", basename(rds_files))
targets <- readLines(targets_file)
targets <- targets[targets != ""] # Remove linhas vazias

sc_objs <- lapply(rds_files, readRDS)

if (seuratObj_path == 'AML_object.rda') {
  load(seuratObj_path)
  seuratObj <- seuratObj[c(VariableFeatures(seuratObj)[1:200], targets[1]),]
} else {
  seuratObj <- readRDS(seuratObj_path)
}

obj_base <- CreateScRank(input = seuratObj, species = species, cell_type = column, target = targets[1])
obj_base@net <- sc_objs
names(obj_base@net) <- cell_types
saveRDS(obj_base, "merged_obj.RDS")

all_ranks <- data.frame()

for (target_sc in targets) {
  message("Processing target: ", target_sc)
  
  tryCatch({
    # A MÁGICA AQUI: Cria o scRank específico para o gene da vez!
    obj <- CreateScRank(input = seuratObj,
                        species = species, 
                        cell_type = column,
                        target = target_sc)
    
    obj@net <- sc_objs
    names(obj@net) <- cell_types
    obj@para$ct.keep = names(obj@net)
    
    obj <- rank_celltype(obj, n.core = 4)
    
    perb_scores <- obj@cell_type_rank$perb_score
    df_long <- data.frame(
      cell_type = cell_types,
      target = target_sc,
      perb_score = as.numeric(perb_scores)
    )
    
    all_ranks <- rbind(all_ranks, df_long)
    message("Finished: ", target_sc)
    
  }, error = function(e) {
    message("Failed for target ", target_sc, ": ", e$message)
  })
}

write.table(
  all_ranks, 
  "perbscore_all_targets.txt", 
  quote = FALSE, 
  row.names = FALSE, 
  col.names = TRUE, 
  sep = "\t"
)
