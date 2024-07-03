######### Plotting correlated CNA-TCs on spatial samples with chromosomes-scale addition #################


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
library(Polychrome)

# Set working directory
setwd("/home/iro/Colocalization_results")
dir()

#### load the correlated CNA-TCs file #########
matching_cna_pairs <- read.table(file = 'matching_cna_pairs.tsv', sep = '\t', header = TRUE)

# read in TACNA-TCs file 
tacna_tcs <- data.frame(read.table("/home/iro/GSEA_results/TACNA_CNA_TCs.tsv", sep = '\t', header = TRUE))

### keep only the chr and tcs columns
tacna_tcs_subset <- tacna_tcs[, c("CHR_Mapping", "Component")]

## rename the component names 
tacna_tcs_subset$Component <- gsub("^([0-9]+)$", "TC\\1", tacna_tcs_subset$Component)


############################################################################
###### Add the chromosomes information in the main dataframe to generate the heatmap 

# Extract the chromosome mappings from tacna_tcs_subset
chr_mappings <- setNames(tacna_tcs_subset$CHR_Mapping, tacna_tcs_subset$Component)

# Use match to find the corresponding chromosomes for tc_1 and tc_2 in matching_cna_pairs
CHR_Mapping_1 <- chr_mappings[match(matching_cna_pairs$tc_1, names(chr_mappings))]
CHR_Mapping_2 <- chr_mappings[match(matching_cna_pairs$tc_2, names(chr_mappings))]

# Add the new columns to matching_cna_pairs using cbind
matching_cna_pairs <- cbind(matching_cna_pairs, CHR_Mapping_1, CHR_Mapping_2)

# Save the dataframe
write.table(matching_cna_pairs, file = "colocalization_results_cna_tcs_with_chromosomes_annotations.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)

# Read the saved dataframe back for further processing
matching_cna_pairs <- read.table(file = 'colocalization_results_cna_tcs_with_chromosomes_annotations.tsv', sep = '\t', header = TRUE)

########## Loop to create heatmaps of the samples for CNA-TCs with correct order of chromosomes #########

# Get unique image IDs
unique_image_ids <- unique(matching_cna_pairs$image_id)

# Create a list to store both initial and wider dataframes for each image id
all_dataframes <- list()

# Extract unique chromosome mappings
unique_chr <- unique(tacna_tcs_subset$CHR_Mapping)

# Define the expected order of chromosomes
chromosome_order <- c(as.character(1:23))

# Ensure unique_chr is in the correct order
unique_chr <- intersect(chromosome_order, unique_chr)

# Define a color palette for chromosomes using the rainbow function
chromosome_colors <- setNames(rainbow(length(unique_chr)), unique_chr)

# Set PDF output for multiple heatmaps
pdf("5colocalization_cna_heatmap.pdf", width = 12, height = 10)

# Loop over unique image IDs
for (image_id in unique_image_ids) {
  # Subset matching_cna_pairs for current image id
  subset_df <- matching_cna_pairs[matching_cna_pairs$image_id == image_id, c("image_id", "tc_1", "tc_2", "coloc_score", "tc_id", "CHR_Mapping_1", "CHR_Mapping_2")]
  
  # Store the initial dataframe in the list with the image id as the list name
  all_dataframes[[paste0("initial_", image_id)]] <- subset_df
  
  # Use dcast to convert subset dataframe into wider format
  wider_df <- dcast(subset_df, tc_1 ~ tc_2, value.var = "coloc_score", fun.aggregate = sum)
  wider_df <- data.frame(wider_df, row.names = 1)
  
  # Remove TCs where all values are zero
  wider_df <- wider_df[rowSums(wider_df != 0) > 0, colSums(wider_df != 0) > 0]
  
  # Extract row and column names for annotations
  row_chromosomes <- subset_df[match(rownames(wider_df), subset_df$tc_1), "CHR_Mapping_1"]
  col_chromosomes <- subset_df[match(colnames(wider_df), subset_df$tc_2), "CHR_Mapping_2"]
  
  # Create annotation dataframes with factor levels set to all unique chromosomes
  row_annotation <- data.frame(Chromosome = factor(row_chromosomes, levels = unique_chr))
  col_annotation <- data.frame(Chromosome = factor(col_chromosomes, levels = unique_chr))
  
  # Set row names for annotations
  rownames(row_annotation) <- rownames(wider_df)
  rownames(col_annotation) <- colnames(wider_df)
  
  # Calculate the maximum absolute value to determine the symmetric color scale
  max_abs <- max(abs(wider_df))
  
  # Define the range of the color scale
  scale_range <- 3
  
  # Define the number of intervals on each side of zero
  num_intervals <- 50
  
  # Calculate the interval size
  interval_size <- scale_range / num_intervals
  
  # Create breaks for the color scale with zero in the center
  breaks <- c(seq(-scale_range, 0, by = interval_size), seq(interval_size, scale_range, by = interval_size))
  
  # Plot heatmap with annotations
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
                clustering_distance_cols = "correlation",
                annotation_row = row_annotation,
                annotation_col = col_annotation,
                annotation_colors = list(Chromosome = chromosome_colors)
  )
  
  # Print the plot
  print(p)
  
  # Store the wider dataframe in the list with the image id as the list name
  all_dataframes[[paste0("wider_", image_id)]] <- wider_df
}

# Close the PDF device
dev.off()

################## loop to create heatmaps of the samples for CNA-TCs using a distinct color palette ##################

# Get unique image IDs
unique_image_ids <- unique(matching_cna_pairs$image_id)

# Create a list to store both initial and wider dataframes for each image id
all_dataframes <- list()

# Extract unique chromosome mappings
unique_chr <- unique(tacna_tcs_subset$CHR_Mapping)

# Define the expected order of chromosomes
chromosome_order <- c(as.character(1:23))

# Ensure unique_chr is in the correct order
unique_chr <- intersect(chromosome_order, unique_chr)

# Generate a 23-color sky palette using Polychrome
sky_palette <- green.armytage.colors(23)[1:length(unique_chr)]
chromosome_colors <- setNames(sky_palette, unique_chr)

# Set PDF output for multiple heatmaps
pdf("7colocalization_cna_heatmap.pdf", width = 12, height = 8)

# Loop over unique image IDs
for (image_id in unique_image_ids) {
  # Subset matching_cna_pairs for current image id
  subset_df <- matching_cna_pairs[matching_cna_pairs$image_id == image_id, c("image_id", "tc_1", "tc_2", "coloc_score", "tc_id", "CHR_Mapping_1", "CHR_Mapping_2")]
  
  # Store the initial dataframe in the list with the image id as the list name
  all_dataframes[[paste0("initial_", image_id)]] <- subset_df
  
  # Use dcast to convert subset dataframe into wider format
  wider_df <- dcast(subset_df, tc_1 ~ tc_2, value.var = "coloc_score", fun.aggregate = sum)
  wider_df <- data.frame(wider_df, row.names = 1)
  
  # Remove TCs where all values are zero
  wider_df <- wider_df[rowSums(wider_df != 0) > 0, colSums(wider_df != 0) > 0]
  
  # Extract row and column names for annotations
  row_chromosomes <- subset_df[match(rownames(wider_df), subset_df$tc_1), "CHR_Mapping_1"]
  col_chromosomes <- subset_df[match(colnames(wider_df), subset_df$tc_2), "CHR_Mapping_2"]
  
  # Create annotation dataframes with factor levels set to all unique chromosomes
  row_annotation <- data.frame(Chromosome = factor(row_chromosomes, levels = unique_chr))
  col_annotation <- data.frame(Chromosome = factor(col_chromosomes, levels = unique_chr))
  
  # Set row names for annotations
  rownames(row_annotation) <- rownames(wider_df)
  rownames(col_annotation) <- colnames(wider_df)
  
  # Calculate the maximum absolute value to determine the symmetric color scale
  max_abs <- max(abs(wider_df))
  
  # Define the range of the color scale
  scale_range <- 3
  
  # Define the number of intervals on each side of zero
  num_intervals <- 50
  
  # Calculate the interval size
  interval_size <- scale_range / num_intervals
  
  # Create breaks for the color scale with zero in the center
  breaks <- c(seq(-scale_range, 0, by = interval_size), seq(interval_size, scale_range, by = interval_size))
  
  # Plot heatmap with annotations
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
                clustering_distance_cols = "correlation",
                annotation_row = row_annotation,
                annotation_col = col_annotation,
                annotation_colors = list(Chromosome = chromosome_colors)
  )
  
  # Print the plot
  print(p)
  
  # Store the wider dataframe in the list with the image id as the list name
  all_dataframes[[paste0("wider_", image_id)]] <- wider_df
}

# Close the PDF device
dev.off()

###### Correct loop to create heatmaps of the samples for CNA-TCs #########

# Get unique image IDs
unique_image_ids <- unique(matching_cna_pairs$image_id)

# Create a list to store both initial and wider dataframes for each image id
all_dataframes <- list()

# Define a color palette for chromosomes
all_chromosomes <- unique(c(matching_cna_pairs$CHR_Mapping_1, matching_cna_pairs$CHR_Mapping_2))
chromosome_colors <- setNames(rainbow(length(all_chromosomes)), all_chromosomes)

# Set PDF output for multiple heatmaps
pdf("6colocalization_cna_heatmap.pdf", width = 12, height = 8)

# Loop over unique image IDs
for (image_id in unique_image_ids) {
  # Subset matching_cna_pairs for current image id
  subset_df <- matching_cna_pairs[matching_cna_pairs$image_id == image_id, c("image_id", "tc_1", "tc_2", "coloc_score", "tc_id", "CHR_Mapping_1", "CHR_Mapping_2")]
  
  # Store the initial dataframe in the list with the image id as the list name
  all_dataframes[[paste0("initial_", image_id)]] <- subset_df
  
  # Use dcast to convert subset dataframe into wider format
  wider_df <- dcast(subset_df, tc_1 ~ tc_2, value.var = "coloc_score", fun.aggregate = sum)
  wider_df <- data.frame(wider_df, row.names = 1)
  
  # Remove TCs where all values are zero
  wider_df <- wider_df[rowSums(wider_df != 0) > 0, colSums(wider_df != 0) > 0]
  
  # Extract row and column names for annotations
  row_chromosomes <- subset_df[match(rownames(wider_df), subset_df$tc_1), "CHR_Mapping_1"]
  col_chromosomes <- subset_df[match(colnames(wider_df), subset_df$tc_2), "CHR_Mapping_2"]
  
  # Create annotation dataframes
  row_annotation <- data.frame(Chromosome = factor(row_chromosomes, levels = all_chromosomes))
  col_annotation <- data.frame(Chromosome = factor(col_chromosomes, levels = all_chromosomes))
  
  # Set row names for annotations
  rownames(row_annotation) <- rownames(wider_df)
  rownames(col_annotation) <- colnames(wider_df)
  
  # Calculate the maximum absolute value to determine the symmetric color scale
  max_abs <- max(abs(wider_df))
  
  # Define the range of the color scale
  scale_range <- 3
  
  # Define the number of intervals on each side of zero
  num_intervals <- 50
  
  # Calculate the interval size
  interval_size <- scale_range / num_intervals
  
  # Create breaks for the color scale with zero in the center
  breaks <- c(seq(-scale_range, 0, by = interval_size), seq(interval_size, scale_range, by = interval_size))
  
  # Plot heatmap with annotations
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
                clustering_distance_cols = "correlation",
                annotation_row = row_annotation,
                annotation_col = col_annotation,
                annotation_colors = list(Chromosome = chromosome_colors)
  )
  
  # Print the plot
  print(p)
  
  # Store the wider dataframe in the list with the image id as the list name
  all_dataframes[[paste0("wider_", image_id)]] <- wider_df
}

# Close the PDF device
dev.off()

