#!/usr/bin/env Rscript
library(Seurat)
library(scRank)
library(dplyr)
library(GENIE3)


args <- commandArgs(trailingOnly = TRUE)

seuratObj <- args[1]
targets <- args[2]
species <- args[3]
column <- args[4]
rds_files <- args[5:length(args)]

cell_types <- sub("_weight.*", "", basename(rds_files))
targets <- readLines(targets)
target <- targets[1]

sc_objs <- lapply(rds_files, readRDS)

if (seuratObj == 'AML_object.rda') {
  load(seuratObj)
  seuratObj <- seuratObj[c(VariableFeatures(seuratObj)[1:200], target),]
} else {
  seuratObj <- readRDS(seuratObj)
}

obj <- CreateScRank(input = seuratObj,
                    species = species, 
                    cell_type = column,
                    target = target)

obj@net <- sc_objs
names(obj@net) <- cell_types

obj@para$ct.keep = names(obj@net)

saveRDS(obj, "merged_obj.RDS")


all_ranks <- data.frame()

for (target_sc in targets) {
  message("Processing target: ", target_sc)
  
  # Set the target
  obj@para$target <- target_sc
  
  # Try running rank_celltype
  tryCatch({
    obj <- rank_celltype(obj, n.core = 4)
    
    # Extract data and convert to long format
    perb_scores <- obj@cell_type_rank$perb_score
    df_long <- data.frame(
      cell_type = cell_types,
      target = target_sc,
      perb_score = as.numeric(perb_scores)
    )
    
    # Append to the main data frame
    all_ranks <- rbind(all_ranks, df_long)
    
    message("Finished: ", target_sc)
    
  }, error = function(e) {
    message("Failed for target ", target_sc, ": ", e$message)
  })
}

# Save results (even if partial)
write.table(
  all_ranks, 
  "perbscore_all_targets.txt", 
  quote = FALSE, 
  row.names = FALSE, 
  col.names = TRUE, 
  sep = "\t"
)


# ===== Final export chunk: GSVA/GSEA/DGE consolidated workbook =====
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("Package 'openxlsx' is required to export Excel. Please install it before running this script.")
}

pick_first_object <- function(candidates) {
  for (nm in candidates) {
    if (exists(nm, inherits = FALSE)) {
      return(get(nm, inherits = FALSE))
    }
  }
  NULL
}

to_df <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x)) return(x)
  if (is.matrix(x)) return(as.data.frame(x, check.names = FALSE))
  as.data.frame(x, check.names = FALSE)
}

parse_comparison_group <- function(source_name) {
  source_name <- as.character(source_name)
  comparison <- source_name
  group <- NA_character_

  # Handle common patterns such as "A_vs_B__Group", "A_vs_B|Group" or "A_vs_B:Group"
  if (grepl("__", source_name, fixed = TRUE)) {
    parts <- strsplit(source_name, "__", fixed = TRUE)[[1]]
    comparison <- parts[1]
    group <- ifelse(length(parts) >= 2, parts[2], NA_character_)
  } else if (grepl("\\|", source_name)) {
    parts <- strsplit(source_name, "\\|", perl = TRUE)[[1]]
    comparison <- parts[1]
    group <- ifelse(length(parts) >= 2, parts[2], NA_character_)
  } else if (grepl(":", source_name, fixed = TRUE)) {
    parts <- strsplit(source_name, ":", fixed = TRUE)[[1]]
    comparison <- parts[1]
    group <- ifelse(length(parts) >= 2, parts[2], NA_character_)
  }

  list(comparison = comparison, group = group)
}

collapse_results_with_labels <- function(x) {
  if (is.null(x)) return(NULL)

  if (is.list(x) && !is.data.frame(x)) {
    nms <- names(x)
    if (is.null(nms) || any(nms == "")) {
      nms <- paste0("comparison_", seq_along(x))
    }

    out_list <- lapply(seq_along(x), function(i) {
      df <- to_df(x[[i]])
      if (is.null(df) || nrow(df) == 0) return(NULL)

      parsed <- parse_comparison_group(nms[i])
      df$comparison <- parsed$comparison
      df$group <- parsed$group
      df
    })
    out_list <- out_list[!vapply(out_list, is.null, logical(1))]
    if (length(out_list) == 0) return(NULL)
    return(dplyr::bind_rows(out_list))
  }

  df <- to_df(x)
  if (is.null(df)) return(NULL)

  if (!"comparison" %in% colnames(df)) {
    df$comparison <- NA_character_
  }
  if (!"group" %in% colnames(df)) {
    df$group <- NA_character_
  }

  df
}

# Candidate objects in case naming differs across scripts
gsva_es <- pick_first_object(c("gsva_es_values", "gsva_es", "gsva_scores", "gsva_score_matrix"))
gsva_diff <- pick_first_object(c("gsva_diff_full", "gsva_diff", "gsva_differential_full", "gsva_differential"))
gsea_all <- pick_first_object(c("gsea_all_results", "gsea_results", "gsea_res", "gsea_full"))
dge_all <- pick_first_object(c("dge_all_results", "dge_results", "deg_results", "dge_full"))

gsva_es_df <- to_df(gsva_es)
gsva_diff_df <- to_df(gsva_diff)
gsea_all_df <- collapse_results_with_labels(gsea_all)
dge_all_df <- collapse_results_with_labels(dge_all)

if (is.null(gsva_es_df)) gsva_es_df <- data.frame(note = "GSVA_es object not found in environment")
if (is.null(gsva_diff_df)) gsva_diff_df <- data.frame(note = "GSVA differential full table not found in environment")
if (is.null(gsea_all_df)) gsea_all_df <- data.frame(note = "GSEA results object not found in environment")
if (is.null(dge_all_df)) dge_all_df <- data.frame(note = "DGE results object not found in environment")

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "GSVA_es_values")
openxlsx::addWorksheet(wb, "GSVA_diff_full")
openxlsx::addWorksheet(wb, "GSEA_all")
openxlsx::addWorksheet(wb, "DGE_all")

openxlsx::writeData(wb, "GSVA_es_values", gsva_es_df)
openxlsx::writeData(wb, "GSVA_diff_full", gsva_diff_df)
openxlsx::writeData(wb, "GSEA_all", gsea_all_df)
openxlsx::writeData(wb, "DGE_all", dge_all_df)

openxlsx::saveWorkbook(wb, "downstream_results.xlsx", overwrite = TRUE)
