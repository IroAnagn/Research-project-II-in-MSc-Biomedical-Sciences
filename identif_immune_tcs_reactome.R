################ Identification of Immune-TCs in Reactome database after GSEA ########################

## Set working directory to the location of GSEA results
setwd("/home/iro/GSEA_results/")
dir()

# Load necessary libraries
library(data.table)
library(dplyr)

# Function to preview a dataframe
pre <- function(df, n_entries=6) {
  rows <- min(nrow(df),n_entries)
  cols <- min(ncol(df),n_entries)
  df[1:rows,1:cols, drop=F]
}


#### load the results from Reactome (c2): 
Reactome <- read.table(file = 'enrichment_matrix_c2.tsv', sep = '\t', header = TRUE)
## set the first column of the data frame as row names:
rownames(Reactome) <- Reactome$X
Reactome <- Reactome[, -1]

# Save the Reactome modified dataframe to a new TSV file
write.table(Reactome, file = "Reactome.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)

# Continue from this step: 
Reactome <- read.table(file = 'Reactome.tsv', sep = '\t', header = TRUE)

########### Bonferroni Correction ####################
Reactome_Bonf <- Reactome

# Function to derive p-values from z-scores
PDerive <- function(Z) {
  P <- pnorm(abs(Z), lower.tail = FALSE) 
  return(P)
}

# Apply PDerive function to each element of the z-score matrix
p_values_Re <- apply(Reactome_Bonf, 2, function(col) PDerive(col))

### histogram of raw p-values 
pdf("hist_raw_p_values_Re.pdf", width = 12, height = 8)
hist(unlist(p_values_Re), breaks=50)
dev.off()

# Apply Bonferroni correction to p-values
p_values_adjusted_Re <- apply(p_values_Re, 2, function(col) p.adjust(col, method = "bonferroni"))

p_values_adjusted_min_Re <- apply(p_values_adjusted_Re, 2, min)

pdf("hist_p_values_React.pdf", width = 12, height = 8)
hist(p_values_adjusted_min_Re, breaks = 40)
dev.off()

table(p_values_adjusted_min_Re < 0.05)

## set the cutoff: p-value < 0.05  
# Find the TCs where the minimum adjusted p-value is less than 0.05
cutoff_TCs_Re <- which(p_values_adjusted_min_Re < 0.05)

# Subset the original data frame to include only the TCs below the cutoff
Reactome_cutoff <- Reactome_Bonf[, cutoff_TCs_Re]

# Check for NAs in the dataframe
any_na <- any(is.na(Reactome_cutoff))
print(any_na)

## hist of the TCs below the cutoff 
# Filter the p_values_adjusted_min_Re vector to include only values below the cutoff
filtered_p_values_Re <- p_values_adjusted_min_Re[p_values_adjusted_min_Re < 0.05]

# Create a histogram of the filtered p-values
pdf("hist_p_values_cutoff_Re.pdf", width = 12, height = 8)
hist(filtered_p_values_Re, breaks = 20, col = "skyblue", 
     main = "Histogram of P-values < 0.05", xlab = "P-value")
dev.off()

################### Ranking of TCs (in Reactome_cutoff) #################
# Apply the rank function to each column/TC and , if ties, keep the lowest p-value 
ranked_Reactome <- apply(Reactome_cutoff, 2, function(col){rank(-abs(col), ties.method = "min")})

# Save the ranked_Reactome dataframe to a TSV file
write.table(ranked_Reactome, file = "ranked_Reactome.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)

# Continue from this step: 
ranked_Reactome <- read.table(file = 'ranked_Reactome.tsv', sep = '\t', header = TRUE)

############# Subsetting only immune GSs using Stefan's gene sets ###############
# load the results from Stefan Reactome: 
Stef_Reactome <- read.table(file = 'reactome_immune_gene_sets_R-HSA-168256_children.tsv', sep = '\t', header = TRUE)

## filter out the FALSE:
Stef_Re_immune <- Stef_Reactome[Stef_Reactome[, 2], ]

## Find the overlapping immune gene sets between Stefan's and my ranked_Reactome data frame ##
# Get the row names of Stef_Re_immune dataframe
Stef_immune_gs_Re <- Stef_Re_immune[, 1] 

# Get the row names of my ranked_Reactome dataframe
ranked_Re_gs <- rownames(ranked_Reactome)

# the gene set names need to be split
split_gs_Re <- strsplit(ranked_Re_gs, " -- ")

# Extract only the 1st part of each split gene set name
split_gs_Re <- lapply(split_gs_Re, `[`, 1)

# Find the overlap now between Stef's and mine 
overlap_gs_Re <- intersect(Stef_immune_gs_Re, split_gs_Re)

############ Choose the immune gene sets and create a matrix with the ranks #######
# Convert the matrix to a data frame
ranked_Re_df <- as.data.frame(ranked_Reactome)

# transform rownames of the df into only the first part of the splitting
rownames(ranked_Re_df) <- split_gs_Re

# keep only the overlapping immune gene sets
subset_ranked_Re_df <- ranked_Re_df[match(overlap_gs_Re, rownames(ranked_Re_df)), , drop = FALSE]
overlap_gs_Re_df <- subset_ranked_Re_df

# Save the overlap_gs_Re_df dataframe as a TSV file
write.table(overlap_gs_Re_df, file = "overlap_gs_Re_df.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)

# Continue from this step: 
overlap_gs_Re_df <- read.table(file = 'overlap_gs_Re_df.tsv', sep = '\t', header = TRUE)

#################### FILTER: Select only the TCs that have a ranking of 1 ######################
# Check if any value in each column is equal to 1
col_rank_1 <- apply(overlap_gs_Re_df[,], 2, function(col) any(col == 1))

# Subset overlap_gs_Re_df based on the columns where at least one TC has a rank of 1
filt_overlap_gs_Re_df <- overlap_gs_Re_df[, col_rank_1]

filt_overlap_gs_Re_df_minimum_ranks = apply(filt_overlap_gs_Re_df,2, min)
summary(filt_overlap_gs_Re_df_minimum_ranks)

write.table(filt_overlap_gs_Re_df, file = "filt_overlap_gs_Re_df.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)

# Continue from this step: 
filt_overlap_gs_Re_df <- read.table(file = 'filt_overlap_gs_Re_df.tsv', sep = '\t', header = TRUE)

# Calculate the column sums
column_sums_Re <- colSums(filt_overlap_gs_Re_df)
column_sums_Re
