##############	creating new mixing matrix #################
### Combine the old CNA-TCs with the new Immune-TCs #### 



# Load required libraries
library(RColorBrewer)
library(Seurat)
library(circlize)
library(data.table)
library(ggplot2)
library(cowplot)

# Set working directory
setwd("/home/iro/Projection_results/")

################ Load the CNA-TCs corrected_mix_mat file ################
cna_subset_corrected_mix_mat <- read.table(file = "cna_subset_corrected_mix_mat.tsv", sep = '\t', header = TRUE)

######## Load the file containing all the TCs ########
all_tcs <- data.frame(read.table("immune_cna_tc_for_projection_spatialsamples.tsv", sep = '\t', header = TRUE)) 

######## Load the correct immune-TCs file ################ 
immune_tcs <- data.frame(read.table("immune_ica_flipped_df.tsv", sep = '\t', header = TRUE)) 

############# Read in all mixing matrices from projection ##############
mixing_matrices <- list.files("/home/iro/Projection_results/", pattern = "mixing_matrix.tsv", full.names = TRUE)
mix_mat_list <- lapply(mixing_matrices, function(file) as.data.frame(fread(file, header = TRUE)))

# Combine the mixing matrices into one dataframe
corrected_mix_mat <- do.call(rbind, mix_mat_list)

# Remove redundant objects from the workspace
rm(mix_mat_list)

# Rename columns in corrected_mix_mat to match the pattern "TC"
colnames(corrected_mix_mat) <- gsub(pattern = "consensus.independent.component.", replacement = "TC", x = colnames(corrected_mix_mat))

#### Merge the CNA- and Immune-corrected mixing matrices ####

# Bind the subsetted corrected_mix_mat with cna_subset_corrected_mix_mat
new_corrected_mixing_matrix <- cbind(cna_subset_corrected_mix_mat, corrected_mix_mat)

# Identify and remove duplicate column names
duplicate_cols <- duplicated(colnames(new_corrected_mixing_matrix))
new_corrected_mixing_matrix <- new_corrected_mixing_matrix[, !duplicate_cols]

# Save the new corrected mixing matrix
write.table(new_corrected_mixing_matrix, file = "new_corrected_mixing_matrix.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)

# Read the new corrected mixing matrix for further processing
new_corrected_mixing_matrix <- read.table(file = "new_corrected_mixing_matrix.tsv", sep = '\t', header = TRUE)

# Read the spatial data 
visium <- readRDS("/home/arkajyotibhattacharya/Projects/spatial_transcriptomics_marco/Data/ALL_VISIUM_FILES_MERGED.rds")
