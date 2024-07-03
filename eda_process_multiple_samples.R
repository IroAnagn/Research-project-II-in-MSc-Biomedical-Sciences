############## Exploratory Data Analysis on all the spatial transcriptomics samples ################

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# Set the working directory to the location of the currently active script (only in RStudio)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Define output directory
OUTPUT_DIR <- getwd()

# Define the path where the samples are stored
samples_path <- "/home/iro/10x_visium_samples"

# Define the pattern to match the sample files
pattern <- "BC|LC|CC|CRC|GB|IDC|LC|M|OC|PC"

# List all files in the directory that match the pattern
files <- list.files(path = samples_path, pattern = pattern, full.names = TRUE)

# Load each dataset individually and assign it to a variable with the sample name
for (file in files) {
  sample_name <- tools::file_path_sans_ext(basename(file))
  assign(sample_name, readRDS(file))
}

# Function to perform EDA on a single sample
process_sample <- function(sample_path) {
  sample_name <- tools::file_path_sans_ext(basename(sample_path))
  sample_output_dir <- file.path(OUTPUT_DIR, sample_name)
  dir.create(sample_output_dir, showWarnings = FALSE)
  setwd(sample_output_dir)
  
  sample_data <- get(sample_name)
  
  # Spatial Feature Plot for nCount_Spatial
  pdf("nCount_spatial.pdf", width = 8, height = 6)
  plot2 <- SpatialFeaturePlot(sample_data, features = "nCount_Spatial") + theme(legend.position = "right")
  print(plot2)
  dev.off()
  
  # Log-normalization and PCA
  mapping_file <- sample_data@images[[sample_name]]@coordinates
  common_samples <- intersect(rownames(mapping_file), colnames(sample_data@assays$Spatial@data))
  mapping_file <- mapping_file[common_samples,]
  raw_counts <- as.matrix(sample_data@assays$Spatial@data[,common_samples])
  
  raw_counts_pca <- raw_counts[rowSums(raw_counts) != 0,]
  raw_counts_pca <- Seurat::LogNormalize(raw_counts_pca)
  
  pca_obj_princomp <- princomp(raw_counts_pca, cor = TRUE, fix_sign = FALSE)
  raw_counts_pc1removed <- pca_obj_princomp$loadings[, -1] %*% t(pca_obj_princomp$scores[, -1])
  
  fwrite(as.data.frame(raw_counts_pc1removed), file.path(getwd(), paste0(sample_name, "_lognorm_pc1removed.tsv")), row.names = TRUE, sep = "\t")
  
  sample_data@assays$Spatial@data <- t(raw_counts_pc1removed)
  saveRDS(sample_data, paste0(sample_name, "_normalized.rds"))
  sample_data <- readRDS(paste0(sample_name, "_normalized.rds"))
  
  # Plot selected genes
  plot_features <- function(sample_data, features, filename) {
    pdf(filename, width = 8, height = 6)
    print(SpatialFeaturePlot(sample_data, features = features))
    dev.off()
  }
  
  plot_features(sample_data, c("ERN2", "FLNA"), paste0(sample_name, "_ERN2_FLNA.pdf"))
  plot_features(sample_data, c("BRCA1", "BRCA2"), paste0(sample_name, "_BRCA1_BRCA2.pdf"))
  plot_features(sample_data, c("HMGA2", "CCNB1"), paste0(sample_name, "_HMGA2_CCNB1.pdf"))
  plot_features(sample_data, c("BUB1", "CDK1"), paste0(sample_name, "_BUB1_CDK1.pdf"))
  
  # Custom plot for BRCA1
  jpeg(filename = file.path(sample_output_dir, "BRCA1_vignette.jpg"), height = 700, width = 1200, quality = 50)
  plot <- SpatialFeaturePlot(sample_data, features = "BRCA1") + 
    theme(legend.text = element_text(size = 0),
          legend.title = element_text(size = 20), legend.key.size = unit(1, "cm"))
  print(plot)
  dev.off()
  
  pdf(paste0(sample_name, "_BRCA1_BRCA1.pdf"), width = 8, height = 6)
  p1 <- SpatialFeaturePlot(sample_data, features = "BRCA1", pt.size.factor = 1.6)
  p2 <- SpatialFeaturePlot(sample_data, features = "BRCA1", alpha = c(0.1, 1))
  print(p1 + p2)
  dev.off()
  
  
  
  ####### Exploratory Data Analysis #################
  counts_matrix <- as.matrix(sample_data@assays$Spatial@counts)
  genes <- rownames(counts_matrix)
  spots <- colnames(counts_matrix)
  
  calculate_summary <- function(values) {
    c(mean(values), median(values), min(values), max(values), quantile(values, 0.25), quantile(values, 0.75))
  }
  
  gene_summary_matrix <- t(apply(counts_matrix, 1, calculate_summary))
  colnames(gene_summary_matrix) <- c("Mean", "Median", "Min", "Max", "1st Quartile", "3rd Quartile")
  
  plot_histograms <- function(summary_df, filename, exclude_zeros = FALSE) {
    pdf(filename, width = 12, height = 8)
    par(mfrow = c(2, 3))
    
    for (stat in colnames(summary_df)) {
      values <- as.numeric(summary_df[, stat])
      if (exclude_zeros) {
        values <- values[values != 0]
      }
      if (length(values) > 0) {
        hist(values, main = paste("Histogram of", stat), xlab = stat, breaks = 40)
      }
    }
    dev.off()
  }
  
  plot_histograms(as.data.frame(gene_summary_matrix), paste0(sample_name, "_gene_summary_histograms.pdf"))
  plot_histograms(as.data.frame(gene_summary_matrix), paste0(sample_name, "_gene_summary_histograms_excl_zeros.pdf"), TRUE)
  
  spots_summary_matrix <- t(apply(counts_matrix, 2, calculate_summary))
  colnames(spots_summary_matrix) <- c("Mean", "Median", "Min", "Max", "1st Quartile", "3rd Quartile")
  
  plot_histograms(as.data.frame(spots_summary_matrix), paste0(sample_name, "_spots_summary_histograms.pdf"))
  plot_histograms(as.data.frame(spots_summary_matrix), paste0(sample_name, "_spots_summary_histograms_excl_zeros.pdf"), TRUE)
}

# Process each sample
for (sample_path in files) {
  process_sample(sample_path)
}
