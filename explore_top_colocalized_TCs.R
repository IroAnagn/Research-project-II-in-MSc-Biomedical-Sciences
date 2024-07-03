############ Explore the biological information of the top colocalized TCs #################



# Load necessary libraries
library(dplyr)  
library(data.table)  

# Read in immune TCs data
immune_tcs <- data.frame(read.table("/home/iro/GSEA_results/immune_ica_flipped_df.tsv", sep = '\t', header = TRUE))

# Read in CNA-TCs data
cna_tcs <- data.frame(read.table("/home/iro/GSEA_results/cna_ica_flipped_df.tsv", sep = '\t', header = TRUE))
rownames(cna_tcs) <- cna_tcs$X  # Set first column as row names
cna_tcs <- cna_tcs[, -1]  # Remove the first column

# Read in combined unique TCs data from immune and CNA datasets
all_tcs <- data.frame(read.table("/home/iro/GSEA_results/immune_cna_tc_for_projection_spatialsamples.tsv", sep = '\t', header = TRUE))

# Load positional data to explore cytobands of CNA-TCs
positional <- data.frame(read.table("/home/iro/GSEA_results/CNA_Pos_df.tsv", sep = '\t', header = TRUE))

# Read in results from GSEA for positional database
enrichment_matrix_c1 <- data.frame(read.table("/home/iro/GSEA_results/enrichment_matrix_c1.tsv", sep = '\t', header = TRUE))

# Load ranked, indicator matrices from GOBP and Reactome databases
gobp <- data.frame(read.table("/home/iro/GSEA_results/filt_overlap_gs_df.tsv", sep = '\t', header = TRUE))
reactome <- data.frame(read.table("/home/iro/GSEA_results/filt_overlap_gs_Re_df.tsv", sep = '\t', header = TRUE))

# Read in results from GSEA for GOBP and Reactome databases
enrichment_matrix_c5 <- data.frame(read.table("/home/iro/GSEA_results/enrichment_matrix_c5.tsv", sep = '\t', header = TRUE))
enrichment_matrix_c2 <- data.frame(read.table("/home/iro/GSEA_results/enrichment_matrix_c2.tsv", sep = '\t', header = TRUE))

# Explore top pathways of immune-TCs
# TC372
rownames(reactome[which(reactome$consensus.independent.component.372 == 1),])

# TC641
rownames(reactome[which(reactome$consensus.independent.component.641 == 1),])

# TC4121
rownames(gobp[which(gobp$consensus.independent.component.4121 == 1),])

# TC8538
rownames(reactome[which(reactome$consensus.independent.component.8538 == 1),])

# TC8628
rownames(reactome[which(reactome$consensus.independent.component.8628 == 1),])

# TC8928
rownames(gobp[which(gobp$consensus.independent.component.8928 == 1),])

# Read in results from ICA to find top genes of top pathways based on gene weights
ica_flipped_tcs <- fread("/home/iro/GSEA_results/ica_flipped_independent_components_consensus.tsv", sep = '\t', header = TRUE)
immune_ica_flipped <- data.frame(read.table("/home/iro/GSEA_results/immune_ica_flipped_df.tsv", sep = '\t', header = TRUE))

# Explore maximum gene weight for EntrezID per TC
# TC372
rownames(immune_ica_flipped)[which.max(immune_ica_flipped$consensus.independent.component.372)]
value <- immune_ica_flipped["3127", "consensus.independent.component.372"]
value

# TC641 
rownames(immune_ica_flipped)[which.max(immune_ica_flipped$consensus.independent.component.641)]
value2 <- immune_ica_flipped["3113", "consensus.independent.component.641"]
value2

# TC4121
rownames(immune_ica_flipped)[which.max(immune_ica_flipped$consensus.independent.component.4121)]
value3 <- immune_ica_flipped["55799", "consensus.independent.component.4121"]
value3

# TC8538 
rownames(immune_ica_flipped)[which.max(immune_ica_flipped$consensus.independent.component.8538)]
value4 <- immune_ica_flipped["2537", "consensus.independent.component.8538"]
value4

# TC8628
rownames(immune_ica_flipped)[which.max(immune_ica_flipped$consensus.independent.component.8628)]
value5 <- immune_ica_flipped["3433", "consensus.independent.component.8628"]
value5

# TC8928 
rownames(immune_ica_flipped)[which.max(immune_ica_flipped$consensus.independent.component.8928)]
value6 <- immune_ica_flipped["440107", "consensus.independent.component.8928"]
value6

# Check gene names of corresponding EntrezID
entrezid_mapping <- fread("/home/iro/GSEA_results/Entrezid_mapping_using_org_Hs_eg_db_17012024.txt", sep = '\t', header = TRUE)

# TC372
TC372 <- entrezid_mapping[entrezid_mapping$mapped_entrez_v1 == 3127, ]

# TC641 
TC641 <- entrezid_mapping[entrezid_mapping$mapped_entrez_v1 == 3113, ]

# TC4121
TC4121 <- entrezid_mapping[entrezid_mapping$mapped_entrez_v1 == 55799, ]

# TC8538 
TC8538 <- entrezid_mapping[entrezid_mapping$mapped_entrez_v1 == 2537, ]

# TC8628
TC8628 <- entrezid_mapping[entrezid_mapping$mapped_entrez_v1 == 3433, ]

# TC8928 
TC8928 <- entrezid_mapping[entrezid_mapping$mapped_entrez_v1 == 440107, ]

# Check top pathways for immune TCs in heatmap example OC2
# TC4962
rownames(reactome[which(reactome$consensus.independent.component.4962 == 1),])

# TC8928
rownames(gobp[which(gobp$consensus.independent.component.8928 == 1),])

# TC2717
rownames(reactome[which(reactome$consensus.independent.component.2717 == 1),])

# TC641 
rownames(reactome[which(reactome$consensus.independent.component.641 == 1),])

# Check highly colocalized pairs of CNA-TCs and immune-TCs across all cancer samples
# Load colocalization scores
cna_immune_tcs_coloc_scores <- data.frame(read.table(file = "/home/iro/Colocalization_results/cnatc1_immunetc2_subset_coloc_scores.tsv", sep = '\t', header = TRUE))

# Filter pairs with coloc_score > 0.7 in all image_ids
above_07_pairs <- cna_immune_tcs_coloc_scores %>%
  group_by(tc_id) %>%
  filter(all(coloc_score > 0.7)) %>%
  distinct(tc_id)

# Filter pairs with coloc_score < -0.4 in all image_ids
below_minus_04_pairs <- cna_immune_tcs_coloc_scores %>%
  group_by(tc_id) %>%
  filter(all(coloc_score < -0.4)) %>%
  distinct(tc_id)

# Find pair with coloc_score > 2 in most cancer samples
above_2_pairs <- cna_immune_tcs_coloc_scores %>%
  filter(coloc_score > 2) %>%
  group_by(tc_id) %>%
  summarise(count = n_distinct(image_id),
            samples = paste(unique(image_id), collapse = ", ")) %>%
  arrange(desc(count)) %>%
  slice(1)

# Define TC pair and cancer samples for correlation scores
tc_pair <- "TC8628_TC8538"
cancer_samples <- c("GB2", "CRC1", "CRC3", "LC2", "M1", "OC1", "BC1", "BC2", "BC3", "BC4", "IDC1")

# Filter for specific TC pair and cancer samples
filtered_data <- cna_immune_tcs_coloc_scores %>%
  filter(tc_id == tc_pair & image_id %in% cancer_samples) %>%
  select(image_id, coloc_score)

# Find pair with coloc_score < -2 in most cancer samples
below_minus_2_pairs <- cna_immune_tcs_coloc_scores %>%
  filter(coloc_score < -2) %>%
  group_by(tc_id) %>%
  summarise(count = n_distinct(image_id),
            samples = paste(unique(image_id), collapse = ", ")) %>%
  arrange(desc(count)) %>%
  slice(1)

# Define TC pair and cancer samples for correlation scores
tc_pair2 <- "TC8914_TC9516"
cancer_samples2 <- c("CRC1", "CRC2", "LC1", "OC1", "OC2", "PC3")

# Filter for specific TC pair and cancer samples
filtered_data <- cna_immune_tcs_coloc_scores %>%
  filter(tc_id == tc_pair2 & image_id %in% cancer_samples2) %>%
  select(image_id, coloc_score)

# Explore specific TCs
# TC9516
rownames(gobp[which(gobp$consensus.independent.component.9516 == 1),])
rownames(immune_ica_flipped)[which.max(immune_ica_flipped$consensus.independent.component.9516)]
value9516 <- immune_ica_flipped["100131187", "consensus.independent.component.9516"]
value9516

TC9516 <- entrezid_mapping[entrezid_mapping$mapped_entrez_v1 == 100131187, ]
