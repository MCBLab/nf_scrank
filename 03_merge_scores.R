#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
target_gene <- args[1]
weight_files <- args[-1]

df_final <- data.frame()

for (f in weight_files) {
  ct <- gsub("_weights\\.rds$", "", basename(f))
  mat <- readRDS(f)
  all_gene_scores <- rowSums(mat)
  gene_ranks <- rank(all_gene_scores) / length(all_gene_scores)
  
  score_val <- ifelse(target_gene %in% names(gene_ranks), gene_ranks[target_gene], 0)
  
  df_final <- rbind(df_final, data.frame(CellType = ct, Score = score_val, Metric = "Rank (%)"))
}

write.csv(df_final, file = paste0(target_gene, "_regulatory_scores.csv"), row.names = FALSE)
message(">>> Scores compilados com sucesso!")
