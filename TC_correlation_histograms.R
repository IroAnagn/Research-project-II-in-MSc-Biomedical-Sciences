############# Histograms on the correlation values of TCs per spatial transcriptomics sample #################



# Set working directory
setwd("/home/iro/Colocalization_results")

# Load necessary libraries
library(reshape2)

############### Load the correlation files ################
# Read in the correlation data for CNA-TCs and immune-TCs
unique_tc_coloc_cna <- data.frame(read.table("unique_tc_coloc_cna.tsv", sep = '\t', header = TRUE))
unique_tc_coloc_immune <- data.frame(read.table("unique_tc_coloc_immune.tsv", sep = '\t', header = TRUE))
tc_coloc_cna_immune <- data.frame(read.table("cnatc1_immunetc2_subset_coloc_scores.tsv", sep = '\t', header = TRUE))

############### Histograms for CNA-TCs ####################
# Exclude the first column (tc_id) when plotting histograms
cna_coloc <- unique_tc_coloc_cna[, -1]

# Calculate the maximum correlation score for each cancer sample
max_corr_cna <- apply(cna_coloc, 2, max)

# Plot bar plot of maximum correlation scores for CNA-TCs across samples
pdf("hist_cna_tcs_coloc.pdf", width = 15, height = 10)
barplot(max_corr_cna, names.arg = colnames(max_corr_cna), 
        main = "Maximum Correlation Scores by Sample", 
        xlab = "Cancer Sample", 
        ylab = "Maximum Correlation Score")
dev.off()

# Plot histograms for individual cancer samples
pdf("hist_cna_samples.pdf", width = 12, height = 8)
histograms <- list()
max_freq <- 80000  # Set the maximum frequency (y-axis limit)
num_breaks <- 30   # Set the number of breaks for the histograms
for (col in colnames(cna_coloc)) {
  hist_data <- hist(cna_coloc[[col]], main = paste("Histogram of", col), xlab = "Correlation values", ylim = c(0, max_freq))
  histograms[[col]] <- hist_data
}
par(mfrow = c(5, 4))  # Adjust the layout according to your number of columns
for (i in 1:length(histograms)) {
  plot(histograms[[i]], main = names(histograms)[i])
}
dev.off()

############### Histograms for Immune-TCs ##################
# Exclude the first column (tc_id) when plotting histograms
immune_coloc <- unique_tc_coloc_immune[, -1]

# Calculate the maximum correlation score for each cancer sample
max_corr_immune <- apply(immune_coloc, 2, max)

# Plot bar plot of maximum correlation scores for immune-TCs across samples
pdf("hist_immune_tcs_coloc.pdf", width = 15, height = 10)
barplot(max_corr_immune, names.arg = colnames(max_corr_immune), 
        main = "Maximum Correlation Scores by Sample", 
        xlab = "Cancer Sample", 
        ylab = "Maximum Correlation Score")
dev.off()

# Plot histograms for individual cancer samples
pdf("hist_immune_samples.pdf", width = 12, height = 8)
histograms <- list()
for (col in colnames(immune_coloc)) {
  hist_data <- hist(immune_coloc[[col]], main = paste("Histogram of", col), xlab = "Correlation values", ylim = c(0, max_freq))
  histograms[[col]] <- hist_data
}
par(mfrow = c(5, 4))  # Adjust the layout according to your number of columns
for (i in 1:length(histograms)) {
  plot(histograms[[i]], main = names(histograms)[i])
}
dev.off()

############## Histograms for Immune-CNA-TCs ###############
# Select relevant columns for immune-CNA-TCs
tc_coloc_cna_immune <- tc_coloc_cna_immune[, c("image_id", "coloc_score", "tc_id")]

# Reshape data using dcast to convert dataframe into wider format
tc_coloc_cna_immune <- dcast(tc_coloc_cna_immune, tc_id ~ image_id, value.var = "coloc_score")

# Exclude the first column (tc_id) when plotting histograms
immune_cna_coloc <- tc_coloc_cna_immune[, -1]

# Calculate the maximum correlation score for each cancer sample
max_corr_immune_cna <- apply(immune_cna_coloc, 2, max)

# Plot bar plot of maximum correlation scores for immune-CNA-TCs across samples
pdf("hist_immune_cna_tcs_coloc.pdf", width = 15, height = 10)
barplot(max_corr_immune_cna, names.arg = colnames(max_corr_immune_cna), 
        main = "Maximum Correlation Scores by Sample", 
        xlab = "Cancer Sample", 
        ylab = "Maximum Correlation Score")
dev.off()

# Plot histograms for individual cancer samples
pdf("hist_immune_cna_samples.pdf", width = 12, height = 8)
histograms <- list()
for (col in colnames(immune_cna_coloc)) {
  hist_data <- hist(immune_cna_coloc[[col]], main = paste("Histogram of", col), xlab = "Correlation values", ylim = c(0, max_freq))
  histograms[[col]] <- hist_data
}
par(mfrow = c(5, 4))  # Adjust the layout according to your number of columns
for (i in 1:length(histograms)) {
  plot(histograms[[i]], main = names(histograms)[i])
}
dev.off()
