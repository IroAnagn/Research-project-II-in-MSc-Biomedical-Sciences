#################### Subsetting the old-wrong corrected mixing matrix ####################

# Set working directory
setwd("/home/iro/Projection_results_wrong")

############## Load the input files ####################
# Load CNA TCs data
cna_tcs <- data.frame(read.table("/home/iro/GSEA_results/cna_ica_flipped_df.tsv", sep = '\t', header = TRUE), row.names = 1) 
# Load immune TCs data
immune_tcs <- data.frame(read.table("/home/iro/GSEA_results/immune_ica_flipped_df.tsv", sep = '\t', header = TRUE)) 
# Load all TCs data
all_tcs <- data.frame(read.table("/home/iro/GSEA_results/immune_cna_tc_for_projection_spatialsamples.tsv", sep = '\t', header = TRUE)) 
# Load old-wrong corrected mixing matrix
corrected_mix_mat <- data.frame(read.table(file = "/home/iro/Projection_results_wrong/corrected_mix_mat.tsv", sep = '\t', header = TRUE))

# Rename columns in cna_tcs to match the pattern "TC"
colnames(cna_tcs) = gsub(pattern = "consensus.independent.component.", replacement = "TC", x = colnames(cna_tcs))

# Subset the corrected_mix_mat to include only CNA-TCs present in both matrices
cna_cols <- intersect(colnames(corrected_mix_mat), colnames(cna_tcs))
cna_subset_corrected_mix_mat <- corrected_mix_mat[, cna_cols]

# Save the subsetted corrected mixing matrix
write.table(cna_subset_corrected_mix_mat, file = "/home/iro/Projection_results_wrong/cna_subset_corrected_mix_mat.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)

# Read the subsetted corrected mixing matrix for further processing
cna_subset_corrected_mix_mat <- read.table(file = "/home/iro/Projection_results_wrong/cna_subset_corrected_mix_mat.tsv", sep = '\t', header = TRUE)
