########### Exploratory Data Analysis on a spatial transcriptomics sample ###
## e.g. OC1

# Load necessary libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# Set working directory
setwd("/home/iro/OC1/")
dir()

# Load the dataset
OC1 <- readRDS("OC1.rds")

# Open a PDF file for exporting nCount_Spatial plot
pdf("nCount_spatial.pdf", width = 8, height = 6)
plot2 <- SpatialFeaturePlot(OC1, features = "nCount_Spatial") + theme(legend.position = "right")
print(plot2)
dev.off()

# Extract mapping coordinates and common samples
image_id <- 'OC1'
mapping_file <- OC1@images[[image_id]]@coordinates
common_samples <- intersect(rownames(mapping_file), colnames(OC1@assays$Spatial@data))
mapping_file <- mapping_file[common_samples,]
raw_counts <- as.matrix(OC1@assays$Spatial@data[,common_samples]) # genes x samples

# LogNormalize the data
raw_counts_pca <- raw_counts[rowSums(raw_counts) != 0,]
raw_counts_pca <- Seurat::LogNormalize(raw_counts_pca)

# Perform PCA and remove PC1
pca_obj_princomp <- princomp(raw_counts_pca, cor = TRUE, fix_sign = FALSE)
raw_counts_pc1removed <- pca_obj_princomp$loadings[,-1] %*% t(pca_obj_princomp$scores[,-1])

# Save the modified data matrix
data.table::fwrite(as.data.frame(raw_counts_pc1removed), file.path(getwd(), paste0(image_id, "_lognorm_pc1removed.tsv")), row.names = TRUE, sep = "\t")

# Update Seurat object with normalized data
OC1_normal <- OC1
OC1_normal@assays$Spatial@data <- t(raw_counts_pc1removed)

# Save and reload normalized Seurat object
saveRDS(OC1_normal, "OC1_normalized.rds")
OC1_normal <- readRDS("OC1_normalized.rds")

# Visualize selected features in the original and normalized datasets
pdf("ERN2_FLNA.pdf", width = 8, height = 6)
SpatialFeaturePlot(OC1, features = c("ERN2", "FLNA"))
SpatialFeaturePlot(OC1_normal, features = c("ERN2", "FLNA"))
dev.off()

pdf("BRCA1_BRCA2.pdf", width = 8, height = 6)
SpatialFeaturePlot(OC1, features = c("BRCA1", "BRCA2"))
SpatialFeaturePlot(OC1_normal, features = c("BRCA1", "BRCA2"))
dev.off()

pdf("HMGA2_CCNB1.pdf", width = 8, height = 6)
SpatialFeaturePlot(OC1, features = c("HMGA2", "CCNB1"))
SpatialFeaturePlot(OC1_normal, features = c("HMGA2", "CCNB1"))
dev.off()

pdf("BUB1_CDK1.pdf", width = 8, height = 6)
SpatialFeaturePlot(OC1, features = c("BUB1", "CDK1"))
SpatialFeaturePlot(OC1_normal, features = c("BUB1", "CDK1"))
dev.off()

# Plot BRCA1 feature with custom theme adjustments
plot <- SpatialFeaturePlot(OC1_normal, features = c("BRCA1")) + 
  theme(legend.text = element_text(size = 0),
        legend.title = element_text(size = 20), 
        legend.key.size = unit(1, "cm"))
jpeg(filename = "/home/arkajyotibhattacharya/bm2/Desktop/Iro_Seurat/OC1/BRCA1_vignette.jpg", height = 700, width = 1200, quality = 50)
print(plot)
dev.off()

# Plot BRCA1 feature with different transparency settings
pdf("BRCA1_BRCA1.pdf", width = 8, height = 6)
p1 <- SpatialFeaturePlot(OC1_normal, features = "BRCA1", pt.size.factor = 1.6)
p2 <- SpatialFeaturePlot(OC1_normal, features = "BRCA1", alpha = c(0.1, 1))
print(p1 + p2)
dev.off()

# Exploratory Data Analysis (EDA)

# Number of spots and genes
num_spots <- ncol(OC1_normal@assays$Spatial@data)
num_genes <- nrow(OC1_normal@assays$Spatial@data)
print(paste("Number of spots:", num_spots))
print(paste("Number of genes:", num_genes))

# Check for missing values in the counts matrix
counts_matrix <- as.matrix(OC1_normal@assays$Spatial@counts)
if (any(is.na(counts_matrix))) {
  print("There are missing values (NA) in the counts matrix.")
} else {
  print("There are no missing values (NA) in the counts matrix.")
}

# Calculate summary statistics for each gene
genes <- rownames(counts_matrix)
gene_summary_matrix <- matrix(NA, nrow = length(genes), ncol = 6)
rownames(gene_summary_matrix) <- genes
colnames(gene_summary_matrix) <- c("Mean", "Median", "Min", "Max", "1st Quartile", "3rd Quartile")

# Function to calculate summary statistics for a gene
calculate_gene_summary <- function(gene_values) {
  c(mean(gene_values), median(gene_values), min(gene_values),
    max(gene_values), quantile(gene_values, 0.25), quantile(gene_values, 0.75))
}

# Apply the function to each gene in counts_matrix
for (gene in genes) {
  gene_values <- counts_matrix[gene, ]
  gene_summary_matrix[gene, ] <- calculate_gene_summary(gene_values)
}
print(gene_summary_matrix)

# Plot histograms for each statistical measure of genes
gene_summary_df <- as.data.frame(gene_summary_matrix)

pdf("gene_summary_histograms.pdf", width = 12, height = 8)
par(mfrow = c(2, 3))
for (stat in colnames(gene_summary_df)) {
  hist(gene_summary_df[, stat], main = paste("Histogram of", stat), xlab = stat)
}
dev.off()

# Plot histograms excluding zero values for each statistical measure of genes
pdf("gene_summary_histograms_excl_zeros.pdf", width = 12, height = 8)
par(mfrow = c(2, 3))
for (stat in colnames(gene_summary_df)) {
  values <- as.numeric(gene_summary_df[, stat])
  if (all(values == 0)) {
    cat("All values are zero for", stat, "\n")
    next
  }
  non_zero_values <- values[values != 0]
  hist(non_zero_values, breaks = 40, main = paste("Histogram of", stat), xlab = stat)
}
dev.off()

# Calculate summary statistics for each spot
spots <- colnames(counts_matrix)
spots_summary_matrix <- matrix(NA, nrow = length(spots), ncol = 6)
rownames(spots_summary_matrix) <- spots
colnames(spots_summary_matrix) <- c("Mean", "Median", "Min", "Max", "1st Quartile", "3rd Quartile")

# Function to calculate summary statistics for a spot
calculate_spots_summary <- function(spot_values) {
  c(mean(spot_values), median(spot_values), min(spot_values),
    max(spot_values), quantile(spot_values, 0.25), quantile(spot_values, 0.75))
}

# Apply the function to each spot in counts_matrix
for (spot in spots) {
  spot_values <- counts_matrix[, spot]
  spots_summary_matrix[spot, ] <- calculate_spots_summary(spot_values)
}
print(spots_summary_matrix)

# Plot histograms for each statistical measure of spots
spot_summary_df <- as.data.frame(spots_summary_matrix)

pdf("spots_summary_histograms.pdf", width = 12, height = 8)
par(mfrow = c(2, 3))
for (stat in colnames(spot_summary_df)) {
  hist(spot_summary_df[, stat], main = paste("Histogram of", stat), xlab = stat)
}
dev.off()

# Plot histograms excluding zero values for each statistical measure of spots
pdf("spots_summary_histograms_excl_zeros.pdf", width = 12, height = 8)
par(mfrow = c(2, 3))
for (stat in colnames(spot_summary_df)) {
  values <- as.numeric(spot_summary_df[, stat])
  if (all(values == 0)) {
    cat("All values are zero for", stat, "\n")
    next
  }
  non_zero_values <- values[values != 0]
  hist(non_zero_values, breaks = 40, main = paste("Histogram of", stat), xlab = stat)
}
dev.off()
