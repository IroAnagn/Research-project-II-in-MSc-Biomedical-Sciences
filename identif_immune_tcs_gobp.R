################### Identification of Immune-TCs in GOBP database after GSEA ########################


# Set working directory to the location of GSEA results
setwd("/home/iro/GSEA_results/")

# Load necessary libraries
library(data.table)
library(dplyr)

# Function to preview a dataframe
pre <- function(df, n_entries = 6) {
  rows <- min(nrow(df), n_entries)
  cols <- min(ncol(df), n_entries)
  df[1:rows, 1:cols, drop = FALSE]
}


# Load the results from GO-BP (c5)
GOBP <- data.frame(fread(file = 'enrichment_matrix_c5.tsv', sep = '\t', header = TRUE), row.names = 1)

############### Bonferroni Correction  ################################
GOBP_Bonf <- GOBP

# Function to derive p-values from z-scores
PDerive <- function(Z) {
  P <- pnorm(abs(Z), lower.tail = FALSE)
  return(P)
}

# Apply PDerive function to each element of the z-score matrix
p_values <- apply(GOBP_Bonf, 2, function(col) PDerive(col))

### Histogram of raw p-values 
pdf("hist_raw_p_values.pdf", width = 12, height = 8)
hist(unlist(p_values), breaks = 50)
dev.off()

# Apply Bonferroni correction to p-values
p_values_adjusted <- apply(p_values, 2, function(col) p.adjust(col, method = "bonferroni"))

# Find the minimum adjusted p-value for each TC
p_values_adjusted_min <- apply(p_values_adjusted, 2, min)

# Create histogram of the minimum adjusted p-values
pdf("hist_p_values.pdf", width = 12, height = 8)
hist(p_values_adjusted_min, breaks = 40)
dev.off()

# Check how many TCs have p value < 0.05
significant_TCs <- table(p_values_adjusted_min < 0.05)

## Set the cutoff: p-value < 0.05  
# Find the TCs where the minimum adjusted p-value is less than 0.05
cutoff_TCs <- which(p_values_adjusted_min < 0.05)

# Subset the original data frame to include only the TCs below the cutoff
GOBP_cutoff <- GOBP[, cutoff_TCs]

# Filter the p_values_adjusted_min vector to include only values below the cutoff
filtered_p_values <- p_values_adjusted_min[p_values_adjusted_min < 0.05]

# Create a histogram of the filtered p-values
pdf("hist_p_values_cutoff.pdf", width = 12, height = 8)
hist(filtered_p_values, breaks = 20, col = "skyblue", 
     main = "Histogram of P-values < 0.05", xlab = "P-value")
dev.off()

################### Ranking of TCs (in GOBP_cutoff) #################
# Apply the rank function to each column/TC and , if ties, keep the lowest p-value 
ranked_GOBP <- apply(GOBP_cutoff, 2, function(col){rank(-abs(col), ties.method = "min")})

# Save the ranked_GOBP dataframe to a TSV file
write.table(ranked_GOBP, file = "ranked_GOBP.tsv", sep = "\t", quote = FALSE)

# Continue from this step: 
ranked_GOBP <- data.frame(fread(file = 'ranked_GOBP.tsv', sep = '\t'), row.names = 1)

############# Subsetting only immune GSs using Stefan's gene sets ###############
# Load the results from Stefan GOBP
Stef_GOBP <- read.table(file = 'gobp_immune_gene_sets_GO0002376_children.tsv', sep = '\t', header = TRUE)

## Filter out the FALSE:
Stef_GOBP_immune <- Stef_GOBP[Stef_GOBP[, 2], ]

## Find the overlapping immune gene sets between Stefan's and my ranked_GOBP data frame ##
# Get the row names of Stef_GOBP_immune dataframe
Stef_immune_gs <- Stef_GOBP_immune[, 1] 

# Get the row names of my ranked_GOBP dataframe
ranked_GOBP_gs <- rownames(ranked_GOBP)

# Split the elements of my character vector (my gene set names in ranked_GOBP_gs)
split_gs <- strsplit(ranked_GOBP_gs, " -- ")

# Extract only the 1st part of each split gene set name
split_gs <- lapply(split_gs, `[`, 1)

# Find the overlap now between Stef's and mine 
overlap_gs <- intersect(Stef_immune_gs, split_gs)

############ Choose the immune gene sets and create a matrix with the ranks #######
# Convert the matrix to a data frame
ranked_GOBP_df <- as.data.frame(ranked_GOBP)

# Transform rownames of the df into only the first part of the splitting
rownames(ranked_GOBP_df) <- split_gs

# Keep only the overlapping immune gene sets
subset_ranked_GOBP_df <- ranked_GOBP_df[match(overlap_gs, rownames(ranked_GOBP_df)), , drop = FALSE]
overlap_gs_df <- subset_ranked_GOBP_df

# Save the overlap_gs_df dataframe as a TSV file
write.table(overlap_gs_df, file = "overlap_gs_df.tsv", sep = "\t", quote = FALSE)

##################### FILTER: Select only the TCs that have a ranking of 1 ######################
# Check if any value in each column is equal to 1
col_rank_1 <- apply(overlap_gs_df[,], 2, function(col) any(col == 1))

# Subset overlap_gs_df based on the columns where at least one TC has a rank of 1
filt_overlap_gs_df <- overlap_gs_df[, col_rank_1]

# Save the indicator matrix to a TSV file
write.table(filt_overlap_gs_df, file = "filt_overlap_gs_df.tsv", sep = "\t", quote = FALSE)

# Calculate the column sums
column_sums_GOBP <- colSums(filt_overlap_gs_df)