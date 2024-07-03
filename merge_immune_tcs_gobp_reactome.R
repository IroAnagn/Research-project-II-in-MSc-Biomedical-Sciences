############################ Immune-TCs Processing Script ################################
# This script processes immune-TCs identified from Reactome and GO Biological Processes (GO-BP).
# It combines and processes these TCs, and generates heatmaps for visual analysis.

# GOBP-TCs stored in: filt_overlap_gs_df
# Reactome-TCs stored in: filt_overlap_gs_Re_df


# Load necessary libraries
library(data.table)
library(RColorBrewer)
library(pheatmap)

# Load input files
filt_overlap_gs_df <- read.table(file = 'path_to_filt_overlap_gs_df', sep = '\t', header = TRUE)
filt_overlap_gs_Re_df <- read.table(file = 'path_to_filt_overlap_gs_Re_df', sep = '\t', header = TRUE)

# Extract column names
TCs_GOBP <- colnames(filt_overlap_gs_df)
TCs_Reactome <- colnames(filt_overlap_gs_Re_df)

# Create a dataframe with unique TCs
TCs_df <- data.frame(TCs = unique(c(TCs_GOBP, TCs_Reactome)))

# Check presence of TCs in GOBP and Reactome databases
TCs_df$GOBP <- TCs_df$TCs %in% TCs_GOBP
TCs_df$Reactome <- TCs_df$TCs %in% TCs_Reactome

# Set the TCs as rownames and remove the TCs column
rownames(TCs_df) <- TCs_df$TCs
TCs_df <- TCs_df[, -1]

# Save the TCs dataframe
write.table(TCs_df, file = "immune_TCs_df.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)

# Load the saved TCs dataframe
TCs_df <- read.table(file = 'immune_TCs_df.tsv', sep = '\t', header = TRUE)

# Load ICA flipped data
ica_flipped_df <- fread(file = 'ica_flipped_independent_components_consensus.tsv', sep = '\t', header = TRUE)
ica_flipped_df <- data.frame(ica_flipped_df, row.names = 1, check.names = TRUE)

# Subset ICA data to keep only immune TCs
rownames_TCs_df <- rownames(TCs_df)
immune_ica_flipped_df <- ica_flipped_df[, rownames_TCs_df, drop = FALSE]

# Save the immune ICA flipped dataframe
write.table(immune_ica_flipped_df, file = "immune_ica_flipped_df.tsv", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

# Load the saved immune ICA flipped dataframe
immune_ica_flipped_df <- read.table(file = 'immune_ica_flipped_df.tsv', sep = '\t', header = TRUE)

################### GOBP Processing ###################
# Subset GOBP_Bonf based on filt_overlap_gs_df
immune_GOBP_df <- GOBP_Bonf[, colnames(filt_overlap_gs_df)]

# Process gene set names in GOBP_Bonf
split_gs_gobp <- strsplit(rownames(immune_GOBP_df), " -- ")
split_gs_gobp <- lapply(split_gs_gobp, `[`, 1)
rownames(immune_GOBP_df) <- split_gs_gobp

# Further subset to keep only matching gene sets
immune_GOBP_df <- immune_GOBP_df[rownames(filt_overlap_gs_df), ]

# Check for NAs
if (any(is.na(immune_GOBP_df))) stop("NA values found in immune_GOBP_df")

# Save the immune GOBP dataframe
write.table(immune_GOBP_df, file = "immune_GOBP_df.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)

# Load the saved immune GOBP dataframe
immune_GOBP_df <- read.table(file = 'immune_GOBP_df.tsv', sep = '\t', header = TRUE)

################### Reactome Processing ###################
# Subset Reactome_Bonf based on filt_overlap_gs_Re_df
immune_Reactome_df <- Reactome_Bonf[, colnames(filt_overlap_gs_Re_df)]

# Process gene set names in Reactome_Bonf
split_gs_reactome <- strsplit(rownames(immune_Reactome_df), " -- ")
split_gs_reactome <- lapply(split_gs_reactome, `[`, 1)
rownames(immune_Reactome_df) <- split_gs_reactome

# Further subset to keep only matching gene sets
immune_Reactome_df <- immune_Reactome_df[rownames(filt_overlap_gs_Re_df), ]

# Check for NAs
if (any(is.na(immune_Reactome_df))) stop("NA values found in immune_Reactome_df")

# Save the immune Reactome dataframe
write.table(immune_Reactome_df, file = "immune_Reactome_df.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)

# Load the saved immune Reactome dataframe
immune_Reactome_df <- read.table(file = 'immune_Reactome_df.tsv', sep = '\t', header = TRUE)

################### Heatmap Generation ###################
# Define color palette
col.pal <- RColorBrewer::brewer.pal(9, "RdYlBu")

# Shorten column names
colnames(immune_Reactome_df) <- gsub("consensus\\.independent\\.component\\.(\\d+)", "TC\\1", colnames(immune_Reactome_df))

# Shorten row names
rownames(immune_Reactome_df) <- gsub("^REACTOME_", "", rownames(immune_Reactome_df))

# Plot heatmap for Reactome data
pdf("heatmap_immune_Reactome_df.pdf", width = 12, height = 8)
pheatmap(immune_Reactome_df, 
         cluster_row = TRUE, 
         cluster_cols = TRUE, 
         color = col.pal, 
         fontsize = 5, 
         fontsize_row = 4, 
         fontsize_col = 4, 
         annotation_legend = FALSE)
dev.off()

# Centering the scale to zero and generating heatmap
max_abs_value <- max(abs(immune_Reactome_df))
num_breaks_each_side <- 20
interval_size <- max_abs_value / num_breaks_each_side
breaks <- c(seq(-max_abs_value, 0, by = interval_size), seq(0, max_abs_value, by = interval_size)[-1])

# Define and reverse the RdYlGn palette
rdylgn_palette <- colorRampPalette(brewer.pal(11, "RdYlGn"))(40)
col.pal <- rev(rdylgn_palette)

pdf("centered_heatmap_immune_Reactome_df.pdf", width = 12, height = 8)
pheatmap(immune_Reactome_df, 
         cluster_row = TRUE, 
         cluster_cols = TRUE, 
         color = col.pal, 
         fontsize = 5, 
         fontsize_row = 4, 
         fontsize_col = 4, 
         annotation_legend = FALSE, 
         breaks = breaks, 
         clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation")
dev.off()

# Process GOBP data similarly
colnames(immune_GOBP_df) <- gsub("consensus\\.independent\\.component\\.(\\d+)", "TC\\1", colnames(immune_GOBP_df))
rownames(immune_GOBP_df) <- gsub("^GOBP_", "", rownames(immune_GOBP_df))

# Centering the scale to zero and generating heatmap for GOBP
max_abs_value <- max(abs(immune_GOBP_df))
interval_size <- max_abs_value / num_breaks_each_side
breaks <- c(seq(-max_abs_value, 0, by = interval_size), seq(0, max_abs_value, by = interval_size)[-1])

pdf("centered_heatmap_immune_GOBP_df.pdf", width = 12, height = 8)
pheatmap(immune_GOBP_df, 
         cluster_row = TRUE, 
         cluster_cols = TRUE, 
         color = col.pal, 
         fontsize = 5, 
         fontsize_row = 4, 
         fontsize_col = 4, 
         annotation_legend = FALSE, 
         breaks = breaks, 
         clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation")
dev.off()
