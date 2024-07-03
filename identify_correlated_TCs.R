######### Plotting correlated TCs on spatial samples #################

# Load necessary libraries
library(Seurat)
library(patchwork)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(data.table)
library(ggplot2)
library(cowplot)
library(reshape2)  
library(pheatmap)
library(tidyr)

# Set working directory
setwd("/home/iro/Colocalization_results")
dir()

### Read in colocalization scores
coloc_scores <- read.table('spatial_kernels_coloc_score4_cutoff2_all.tsv', sep = '\t', header = TRUE)

# Read in immune and CNA TCs
immune_tcs <- read.table("/home/iro/GSEA_results/immune_ica_flipped_df.tsv", sep = '\t', header = TRUE)
cna_tcs <- read.table("/home/iro/GSEA_results/cna_ica_flipped_df.tsv", sep = '\t', header = TRUE)

# Set the first column as row names for cna_tcs and remove it
rownames(cna_tcs) <- cna_tcs$X
cna_tcs <- cna_tcs[, -1]

# Identify common TCs
common_tcs <- intersect(colnames(immune_tcs), colnames(cna_tcs))

# Read all TCs combined
all_tcs <- read.table("/home/iro/GSEA_results/immune_cna_tc_for_projection_spatialsamples.tsv", sep = '\t', header = TRUE)

# Rename columns to standard format
immune_comp_numbers <- gsub(".*\\.(\\d+)$", "\\1", colnames(immune_tcs))
cna_comp_numbers <- gsub(".*\\.(\\d+)$", "\\1", colnames(cna_tcs))

colnames(immune_tcs) <- paste0("TC", immune_comp_numbers)
colnames(cna_tcs) <- paste0("TC", cna_comp_numbers)

# Subset matching pairs in colocalization scores
matching_cna_pairs <- coloc_scores[
  coloc_scores$tc_1 %in% colnames(cna_tcs) & coloc_scores$tc_2 %in% colnames(cna_tcs),
]
matching_immune_pairs <- coloc_scores[
  coloc_scores$tc_1 %in% colnames(immune_tcs) & coloc_scores$tc_2 %in% colnames(immune_tcs),
]

# Save matching pairs
write.table(matching_cna_pairs, file = "matching_cna_pairs.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(matching_immune_pairs, file = "matching_immune_pairs.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)

# Load saved matching pairs
matching_cna_pairs <- read.table('matching_cna_pairs.tsv', sep = '\t', header = TRUE)
matching_immune_pairs <- read.table('matching_immune_pairs.tsv', sep = '\t', header = TRUE)

# Generate heatmaps for CNA TCs
unique_image_ids <- unique(matching_cna_pairs$image_id)
all_dataframes <- list()

pdf("4colocalization_cna_heatmap.pdf", width = 12, height = 8)
for (image_id in unique_image_ids) {
  subset_df <- matching_cna_pairs[matching_cna_pairs$image_id == image_id, c("image_id", "tc_1", "tc_2", "coloc_score", "tc_id")]
  all_dataframes[[paste0("initial_", image_id)]] <- subset_df
  
  wider_df <- dcast(subset_df, tc_1 ~ tc_2, value.var = "coloc_score", fun.aggregate = sum)
  wider_df <- data.frame(wider_df, row.names = 1)
  wider_df <- wider_df[rowSums(wider_df != 0) > 0, colSums(wider_df != 0) > 0]
  
  max_abs <- max(abs(wider_df))
  scale_range <- 3
  num_intervals <- 50
  interval_size <- scale_range / num_intervals
  breaks <- c(seq(-scale_range, 0, by = interval_size), seq(interval_size, scale_range, by = interval_size))
  
  p <- pheatmap(wider_df, 
                cluster_row = TRUE,
                cluster_cols = TRUE,
                fontsize = 3,
                fontsize_row = 3, 
                fontsize_col = 3,
                annotation_legend = TRUE,
                main = image_id, 
                breaks = breaks,
                clustering_distance_rows = "correlation",  
                clustering_distance_cols = "correlation"
  )
  
  print(p)
  all_dataframes[[paste0("wider_", image_id)]] <- wider_df
}
dev.off()

# Generate heatmaps for Immune TCs
unique_image_ids_immune <- unique(matching_immune_pairs$image_id)
all_dataframes_immune <- list()

pdf("4colocalization_immune_heatmap.pdf", width = 12, height = 10)
for (image_id in unique_image_ids_immune) {
  subset_df_immune <- matching_immune_pairs[matching_immune_pairs$image_id == image_id, c("image_id", "tc_1", "tc_2", "coloc_score", "tc_id")]
  all_dataframes_immune[[paste0("initial_", image_id)]] <- subset_df_immune
  
  wider_df <- dcast(subset_df_immune, tc_1 ~ tc_2, value.var = "coloc_score", fun.aggregate = sum)
  wider_df <- data.frame(wider_df, row.names = 1)
  wider_df <- wider_df[rowSums(wider_df != 0) > 0, colSums(wider_df != 0) > 0]
  
  max_abs <- max(abs(wider_df))
  scale_range <- 3
  num_intervals <- 50
  interval_size <- scale_range / num_intervals
  breaks <- c(seq(-scale_range, 0, by = interval_size), seq(interval_size, scale_range, by = interval_size))
  
  p <- pheatmap(wider_df, 
                cluster_row = TRUE,
                cluster_cols = TRUE,
                fontsize = 3,
                fontsize_row = 3, 
                fontsize_col = 3,
                annotation_legend = TRUE,
                main = image_id, 
                breaks = breaks, 
                clustering_distance_rows = "correlation",  
                clustering_distance_cols = "correlation"
  )
  
  print(p)
  all_dataframes_immune[[paste0("wider_", image_id)]] <- wider_df
}
dev.off()

# Identify top TC pairs
subset2_cna <- matching_cna_pairs[, c("image_id", "coloc_score", "tc_id")]
subset2_immune <- matching_immune_pairs[, c("image_id", "coloc_score", "tc_id")] 

tc_coloc_cna <- dcast(subset2_cna, tc_id ~ image_id, value.var = "coloc_score")
tc_coloc_immune <- dcast(subset2_immune, tc_id ~ image_id, value.var = "coloc_score")

split_tc_ids <- strsplit(tc_coloc_cna$tc_id, "_")
split_tc_ids_immune <- strsplit(tc_coloc_immune$tc_id, "_")

unique_tc_pairs <- tc_coloc_cna[sapply(split_tc_ids, function(x) length(unique(x)) > 1), ]
unique_tc_pairs_immune <- tc_coloc_immune[sapply(split_tc_ids_immune, function(x) length(unique(x)) > 1), ]

write.table(unique_tc_pairs, file = "unique_tc_coloc_cna.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(unique_tc_pairs_immune, file = "unique_tc_coloc_immune.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)

unique_tc_coloc_cna <- read.table('unique_tc_coloc_cna.tsv', sep = '\t', header = TRUE)
unique_tc_coloc_immune <- read.table('unique_tc_coloc_immune.tsv', sep = '\t', header = TRUE)

max_cna_corr <- max(unique_tc_coloc_cna[, 2:22], na.rm = TRUE)
min_cna_corr <- min(unique_tc_coloc_cna[, 2:22], na.rm = TRUE)
max_immune_corr <- max(unique_tc_coloc_immune[, 2:22], na.rm = TRUE)
min_immune_corr <- min(unique_tc_coloc_immune[, 2:22], na.rm = TRUE)

unique_tc_coloc_cna <- unique_tc_coloc_cna %>%
  separate(tc_id, into = c("tc_1", "tc_2"), sep = "_")  
unique_tc_coloc_immune <- unique_tc_coloc_immune %>%
  separate(tc_id, into = c("tc_1", "tc_2"), sep = "_")  

unique_tc_coloc_cna_tc1 <- unique_tc_coloc_cna$tc_1
unique_tc_coloc_immune_tc2 <- unique_tc_coloc_immune$tc_2

subset_coloc_scores <- coloc_scores %>%
  filter(tc_1 %in% unique_tc_coloc_cna_tc1 & tc_2 %in% unique_tc_coloc_immune_tc2) %>%
  filter(tc_1 != tc_2)

write.table(subset_coloc_scores, file = "cnatc1_immunetc2_subset_coloc_scores.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
subset_coloc_scores <- read.table('cnatc1_immunetc2_subset_coloc_scores.tsv', sep = '\t', header = TRUE)

max_cna_immune_corr <- max(subset_coloc_scores$coloc_score)
min_cna_immune_corr <- min(subset_coloc_scores$coloc_score)
