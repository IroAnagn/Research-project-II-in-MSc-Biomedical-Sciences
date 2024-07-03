############# Identification of CNA-TCs ########################

# Load results from Gene Set Enrichment Analysis (GSEA) for Positional Data (c1): 

Positional <- read.table(file = 'enrichment_matrix_c1.tsv', sep = '\t', header = TRUE)
## set the first column of the data frame as row names:
rownames(Positional) <- Positional$X
Positional <- Positional[, -1]

# Save the Positional modified dataframe to a new TSV file
write.table(Positional, file = "Positional_modif.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)

# Continue from this step: 
Positional <- read.table(file = 'Positional_modif.tsv', sep = '\t', header = TRUE)


################### Bonferroni Correction #####################
Positional_Bonf <- Positional

# Function to derive p-values from z-scores
PDerive <- function(Z) {
  P <- pnorm(abs(Z), lower.tail = FALSE) 
  return(P)
}

# Apply PDerive function to each element of the z-score matrix
p_values_Pos <- apply(Positional_Bonf, 2, function(col) PDerive(col))
#  otherwise:  p_values <- apply(Positional_Bonf, 2, PDerive)


### histogram of raw p-values 
pdf("hist_raw_p_values_Pos.pdf", width = 12, height = 8)
hist(unlist(p_values_Pos), breaks=50)
dev.off()


# Apply Bonferroni correction to p-values
p_values_adjusted_Pos <- apply(p_values_Pos, 2, function(col) p.adjust(col, method = "bonferroni"))
p_values_adjusted_min_Pos = apply(p_values_adjusted_Pos, 2, min)

pdf("hist_p_values_Pos.pdf", width = 12, height = 8)
hist(p_values_adjusted_min_Pos, breaks = 40)
dev.off()

table(p_values_adjusted_min_Pos < 0.05)

## set the cutoff: p-value < 0.05  
# Find the TCs where the minimum adjusted p-value is less than 0.05
cutoff_TCs_Pos <- which(p_values_adjusted_min_Pos < 0.05)

# Subset the original data frame to include only the TCs below the cutoff
Positional_cutoff <- p_values_adjusted_Pos[, cutoff_TCs_Pos]

## hist of the TCs below the cutoff 
# Filter the p_values_adjusted_min_Pos vector to include only values below the cutoff
filtered_p_values_Pos <- p_values_adjusted_min_Pos[p_values_adjusted_min_Pos < 0.05]

# Create a histogram of the filtered p-values
pdf("hist_p_values_cutoff_Pos.pdf", width = 12, height = 8)
hist(filtered_p_values_Pos, breaks = 20, col = "skyblue", 
     main = "Histogram of P-values < 0.05", xlab = "P-value")

dev.off()

########## if p-value < 0.05 convert to 1, and if p-value > 0.05 convert to 0 ##############
########################## Create the indicator matrix #################
Positional_indicator <- ifelse(Positional_cutoff < 0.05, 1, 0)

# Column sum 
pdf("hist_colSums_Pos.pdf", width = 12, height = 8)
hist(colSums(Positional_indicator))
dev.off()

length(which(colSums(Positional_indicator) >1))
sort(colSums(Positional_indicator))

Positional_indicator <- Positional_indicator[, which(colSums(Positional_indicator) >1)]

# Convert the matrix to a dataframe
Positional_indicator_df <- as.data.frame(Positional_indicator)

## Save the indicator df in a tsv file 
write.table(Positional_indicator_df, file = "Positional_indicator.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
# Continue from this step: 
Positional_indicator_df <- read.table(file = 'Positional_indicator.tsv', sep = '\t', header = TRUE)


############################ Heatmap of the CNA-TCs ##########################################################
## CNA-TCs with z-scores before Bonferroni Correction & after GSEA: 359 CNA-TCs 

# Keep only the 359 CNA-TCs in Positional_Bonf and the 269 cytobands 
CNA_Pos <- intersect(colnames(Positional_indicator_df), colnames(Positional_Bonf))
CNA_Pos_df <- Positional_Bonf[rownames(Positional_Bonf), CNA_Pos]

# save it 
write.table(CNA_Pos_df, file = "CNA_Pos_df.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
# Continue from this step: 
CNA_Pos_df <- read.table(file = 'CNA_Pos_df.tsv', sep = '\t', header = TRUE)

(col.pal <- RColorBrewer::brewer.pal(9, "RdYlBu"))
pacman::p_load(pheatmap)

CNA_Pos_df2 <- CNA_Pos_df
# split the names of the cytobands of CNA_Pos_df
split_CNA <- strsplit(rownames(CNA_Pos_df2), " -- ")
# Extract the first part 
first_part_Pos <- sapply(split_CNA, "[", 1)
# Assign the first part as new row names to the CNA_Pos_df
rownames(CNA_Pos_df2) <- first_part_Pos


# Transform column names
colnames(CNA_Pos_df2) <- gsub("consensus\\.independent\\.component\\.(\\d+)", "comp.\\1", 
                              colnames(CNA_Pos_df2))
# save it 
write.table(CNA_Pos_df2, file = "CNA_Pos_df2.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
# Continue from this step: 
CNA_Pos_df2 <- read.table(file = 'CNA_Pos_df2.tsv', sep = '\t', header = TRUE)


pdf("heatmap_CNA_Pos_df2.pdf", width = 12, height = 8)

pheatmap::pheatmap((CNA_Pos_df2), 
                   cluster_row = T,
                   cluster_cols = T,
                   color = col.pal, 
                   fontsize = 5,
                   fontsize_row = 4, 
                   fontsize_col = 4,
                   annotation_legend = FALSE,  # Do not show separate legend for annotations
                   # labels_row = rep(" ", nrow(subset_counts_matrix_1)),
                   # labels_col = rep(" ", ncol(subset_counts_matrix_1)),
                   gaps_col=50)

# Close the PDF device
dev.off()


####### load the TACNA results from Arko #################
TACNA <- read.table(file = '/home/iro/Genomic mapping/Genelevel_DEGR_TCGA_100_EV_All_Components_Amplified_or_Deleted_details_FDR_0.05_CL_0.5_state_deciding_cutoff_0.85_.txt', sep = "\t", header = TRUE)

## TACNA TCs in col: Component 
# Get the unique TCs along with their frequency count
unique_TCs <- table(TACNA$Component)
# Convert the result to a dataframe for better readability
unique_TCs_df <- as.data.frame(unique_TCs)
TACNA_subset <- unique_TCs_df$Freq > 10
sum(TACNA_subset)

TACNA_sub <- unique_TCs_df[TACNA_subset, ]
TACNA_sub$Var2 <- paste( "comp.", TACNA_sub$Var1, sep = "", collapse = NULL)

# Subset the original TACNA file for the CNA-TCs: 
# Find the common TCs between TACNA and TACNA_sub
common_TCs <- intersect(TACNA$Component, TACNA_sub$Var1)

# Subset TACNA to keep only rows with common TCs
TACNA_CNA_TCs <- TACNA[TACNA$Component %in% common_TCs, ]
unique_TACNA_CNA_TCs <- table(TACNA_CNA_TCs$Component)

# save it 
write.table(TACNA_CNA_TCs, file = "TACNA_CNA_TCs.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
# Continue from this step: 
TACNA_CNA_TCs <- read.table(file = 'TACNA_CNA_TCs.tsv', sep = '\t', header = TRUE)

################### Subset the ICA_flipped file for the 476 CNA-TCs from TACNA ########
# load it 
ica_flipped <- fread(file = 'ica_flipped_independent_components_consensus.tsv', sep = '\t', header = TRUE)

# Transform the names of components so they match 
TACNA_sub$Var3 <- paste( "consensus.independent.component.", TACNA_sub$Var1, sep = "", collapse = NULL)

ica_flipped_df <- data.frame(ica_flipped, row.names = 1, check.names=T)
pre(ica_flipped_df)

## keeping only 476 CNA-TCs
# Extract the rownames from TCs_df
CNA_TCs <- TACNA_sub$Var3

# Subset the ica_flipped_df based on matching column names with TCs_df
cna_ica_flipped_df <- ica_flipped_df[, CNA_TCs, drop = FALSE]

# Save it
write.table(cna_ica_flipped_df, file = "cna_ica_flipped_df.tsv", sep = "\t", row.names = TRUE, col.names = NA, quote = F)
# Continue from this step: 
cna_ica_flipped_df <- read.table(file = 'cna_ica_flipped_df.tsv', sep = '\t', header = TRUE)

################### Heatmap of TACNA-CNA-TCs ############################
#### create a dataframe out of TACNA_sub and Positional #########

# first, make smaller row names in Positional 
split_rows_pos <- strsplit(rownames(Positional), " -- ")
head(split_rows_pos)

# Extract only the 1st part of each split gene set name
split_rows_pos <- lapply(split_rows_pos, `[`, 1)
head(split_rows_pos)

# make smaller col names in Positional
split_cols_pos <- gsub("consensus\\.independent\\.component\\.(\\d+)", "comp.\\1", colnames(Positional))
head(split_cols_pos)

Positional2 <- Positional

##### new Positional dataframe with smaller names 
rownames(Positional2) <- split_rows_pos
colnames(Positional2) <- split_cols_pos

# Save Positional2
write.table(Positional2, file = "Positional_shorter_names.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
# Continue from this step: 
Positional2 <- read.table(file = 'Positional_shorter_names.tsv', sep = '\t', header = TRUE)

##### now create a dataframe out of TACNA_sub and Positional2
tcs_pos_and_tacna <- intersect(colnames(Positional2), TACNA_sub$Var2)
tacna_df_zscores <- Positional2[, tcs_pos_and_tacna]

# Save it 
write.table(tacna_df_zscores, file = "tacna_df_zscores.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
# Continue from this step: 
tacna_df_zscores <- read.table(file = 'tacna_df_zscores.tsv', sep = '\t', header = TRUE)

library(pheatmap)
# Define a custom color palette with white as minimum and red as maximum
my_palette <- colorRampPalette(c("white", "red"))(100)

# Find the maximum absolute z-score value
max_abs_z <- max(abs(tacna_df_zscores))
max_abs_z
# heatmap
pdf("heatmap2_TACNA_CNA_TCs.pdf", width = 12, height = 8)

pheatmap::pheatmap((tacna_df_zscores), 
                   cluster_row = T,
                   cluster_cols = T,
                   color = col.pal, 
                   fontsize = 5,
                   fontsize_row = 4, 
                   fontsize_col = 4,
                   annotation_legend = FALSE,  # Do not show separate legend for annotations
                   # labels_row = rep(" ", nrow(subset_counts_matrix_1)),
                   # labels_col = rep(" ", ncol(subset_counts_matrix_1)),
                   gaps_col=50)

# Close the PDF device
dev.off()

#### repeat it but with distance: 1-cor instead of Euclidean 
library(RColorBrewer)
# Define the number of colors in the palette
num_colors <- 40

# Get the "RdYlGn" color palette
rdylgn_palette <- colorRampPalette(brewer.pal(11, "RdYlGn"))(num_colors)
# Reverse the RdYlGn palette
rdylgn_palette_inverted <- rev(rdylgn_palette)

col.pal <- rdylgn_palette_inverted

# Find the maximum absolute z-score value
max_abs_z <- max(abs(tacna_df_zscores))
max_abs_z

# Determine the number of breaks on each side of zero
num_breaks_each_side <- 20

# Calculate the interval size
interval_size <- max_abs_value / num_breaks_each_side

# Generate breaks symmetrically around zero
breaks <- c(seq(-max_abs_value, 0, by = interval_size),
            seq(0, max_abs_value, by = interval_size)[-1])

# Plot the heatmap using the final palette and breaks
pdf("3heatmap_TACNA_CNA_TCs.pdf", width = 12, height = 8)

pheatmap::pheatmap(tacna_df_zscores, 
                   cluster_row = TRUE,
                   cluster_cols = TRUE,
                   color = col.pal, 
                   fontsize = 5,
                   fontsize_row = 4, 
                   fontsize_col = 4,
                   annotation_legend = FALSE,
                   breaks = breaks,
                   gaps_col = 50, 
                   clustering_distance_rows = "correlation",  
                   clustering_distance_cols = "correlation"
)

# Close the PDF device
dev.off()



######### Create	1 matrix with all unique TCs of all 3 datasets (immune- & CNA- TCs)
### combine the cna_ica_flipped_df and immune_ica_flipped_df but only the unique TCs 

# read in the immune_ica_flipped_df file 
immune_ica_flipped_df <-  data.frame(read.table("/home/iro/GSEA_results/immune_ica_flipped_df.tsv", sep = '\t', header = TRUE))

# Extract TC columns from both data frames
cna_tc_cols <- colnames(cna_ica_flipped_df)[-1]  
immune_tc_cols <- colnames(immune_ica_flipped_df)  

common <- intersect(cna_tc_cols, immune_tc_cols)
# there are 15 common TCs between CNAs and Immune 

# Extract the first column as row names
rownames(cna_ica_flipped_df) <- cna_ica_flipped_df$X

# Remove the first column from the data frame
cna_ica_flipped_df <- cna_ica_flipped_df[, -1]

df_merged <- cbind(cna_ica_flipped_df, immune_ica_flipped_df)
df_merged$consensus.independent.component.40

# Identify and remove duplicate column names
duplicate_cols <- duplicated(colnames(df_merged))
# Extract duplicate column names
duplicate_names <- colnames(df_merged)[duplicate_cols]
duplicate_names

new_df_merged <- df_merged[, !duplicate_cols]

write.table(new_df_merged, file = "immune_cna_tc_for_projection_spatialsamples.tsv", sep = "\t", quote = F, row.names = TRUE, col.names = TRUE)

# Continue from this step: 
new_df_merged <- read.table(file = 'immune_cna_tc_for_projection_spatialsamples.tsv', sep = '\t', header = TRUE)

