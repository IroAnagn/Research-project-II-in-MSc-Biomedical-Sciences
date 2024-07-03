Readme file for R scripts generated during Research Project II.

“Investigating the relation between copy number alterations and immune transcriptional footprints in the context of spatial transcriptomics in cancer”

(WMBM902-30)


Intern: Argyro Anagnostopoulou (S5340411)
Daily supervisor(s): dr. A. (Arkajyoti) Bhattacharya, S. (Stefan) Loipfinger
Supervisor: Prof. dr. R.S.N. (Rudolf) Fehrmann
First examiner: Prof. dr. M. (Marco) de Bruyn
MSc Biomedical Sciences: Biology of Cancer and Immune System (2022-2024)
Translational Bioinformatics Group, Department of Medical Oncology, University Medical Center Groningen (UMCG), the Netherlands
Start date: 15/01/2024 - End date: 16/06/2024



#####################################################################################

Data generated during this research project are available in the website https://sites.google.com/student.rug.nl/research-project2-s5340411. 


How to run:
The scripts are ordered in the way the analysis was executed: 
1. eda_spatial_sample.R 
2. eda_process_multiple_samples.R
3. identif_immune_tcs_reactome.R
4. identif_immune_tcs_gobp.R
5. identif_cna_tcs.R 
6. merge_immune_tcs_gobp_reactome.R
7. genomic_mapping.R
8. plotting_tcs_on_spatial_samples.R
9. plotting_CNA_TCs_on_spatial_samples.R
10. check_colocalization_ks_slurm_part1.R
11. check_colocalization_ks_slurm_part2.R
12. check_colocalization_ks_slurm_part3.R
13. identify_correlated_TCs.R
14. plotting_top_colocalized_tcs_on_spatial_samples.R
15. TC_correlation_histograms.R
16. subset_old_corrected_mixing_matrix.R
17. creating_new_mixing_matrix.R
18. cna_tcs_heatmap_chr_annotations.R
19. explore_top_colocalized_TCs


###################################################################################

Script: eda_spatial_sample.R 

Description: 
This R script is for performing exploratory data analysis (EDA) on spatial transcriptomic samples using the Seurat package. The script includes steps for data normalization, principal component analysis (PCA), and visualization of spatial feature plots. It is tailored for the sample OC1 but can be adapted for other samples.

Input Files
A Seurat object file (e.g., sample.rds) containing the spatial transcriptomic data.

Output Files
PDF files with spatial feature plots for various genes.
TSV file with log-normalized counts after removing the first principal component.
RDS file with the normalized data.
PDF files with histograms of summary statistics (mean, median, min, max, 1st quartile, 3rd quartile) for each gene and spot, including versions that exclude zero values.
JPEG image with a customized plot for a specific feature.


#######################################################################################

Script: eda_process_multiple_samples.R

Description: 
This R script is for performing EDA on multiple Visium spatial transcriptomics samples using the Seurat package. The script reads in multiple datasets, performs normalization and PCA, and generates various plots and summary statistics for each sample.

Input Files
Spatial Transcriptomics Sample Files: .rds files containing the spatial transcriptomics data for different samples. These files should be stored in the directory specified by the samples_path variable in the script.

Output Files
For each sample, the script generates the following outputs in a separate directory named after the sample:
Spatial Feature Plots: PDFs and images of spatial feature plots for various genes.
Normalized Data: Log-normalized data with the first principal component removed (.tsv and .rds files).
Summary Statistics and Histograms: PDFs of histograms for gene and spot expression statistics.



########################################################################################

Script: identif_immune_tcs_reactome.R

Description: 
This R script is for identifying TCs associated with immune processes (immune-TCs) in the Reactome database using results obtained through Gene Set Enrichment Analysis (GSEA). The script includes steps for Bonferroni correction, ranking of TCs, and selection of immune-related gene sets.

Data Processing Steps:
* Loading Data: Read the enrichment matrix file from the Gene Ontology Biological Processes (GOBP) database.
* Bonferroni Correction: Apply Bonferroni correction to p-values derived from z-scores.
* Histogram Analysis: Visualize raw and adjusted p-values to assess significance.
* Selection of Significant TCs: Identify TCs with adjusted p-values below the significance threshold (p < 0.05).
* Ranking of TCs: Rank TCs based on their absolute z-scores.
* Filtering Immune Gene Sets: Use immune gene sets from external sources to select TCs associated with immune-related processes.
* Output of Selected TCs: Save the final selected TCs, ranked 1, for further analysis.

Input Files
Tab-separated file containing Reactome pathway enrichment results from GSEA.

Output Files
Tab-separated file containing the filtered and ranked TCs associated with immune processes. 


#######################################################################################

Script: identif_immune_tcs_gobp.R

Description: 
This R script is for identifying TCs associated with immune-related gene sets in the GOBP database using results from GSEA.

Data Processing Steps:
* Loading Data: Read the enrichment matrix file from the GOBP database.
* Bonferroni Correction: Apply Bonferroni correction to p-values derived from z-scores.
* Histogram Analysis: Visualize raw and adjusted p-values to assess significance.
* Selection of Significant TCs: Identify TCs with adjusted p-values below the significance threshold (p < 0.05).
* Ranking of TCs: Rank TCs based on their absolute z-scores.
* Filtering Immune Gene Sets: Use immune gene sets from external sources to select TCs associated with immune-related processes.
* Output of Selected TCs: Save the final selected TCs, ranked 1, for further analysis.

Input Files:
Tab-separated values file containing enrichment matrix from the GOBP database.

Output Files:
Tab-separated values file containing the selected TCs associated with immune-related gene sets, ranked 1.


#######################################################################################

Script: identif_cna_tcs.R 

Description:
This R script is for identifying CNA-TCs using results from Gene Set Enrichment Analysis (GSEA) and the results from the Transcriptional Adaptation to CNA profiling (TACNA) methodology and performing downstream analysis including the generation of heatmaps.

Data Processing Steps:
* Load GSEA Results
* Bonferroni Correction
* Identify Significant TCs
* Generate Heatmaps
* Subset TACNA and ICA Results
* Merge Dataframes

Input Files:
* GSEA Results: Enrichment matrix in TSV format.
* TACNA Results: TCs and associated information in TSV format.
* ICA Flipped Results: TCs and gene expression data in TSV format.

Output Files:
* Positional_modif.tsv: Modified positional dataframe.
* hist_raw_p_values_Pos.pdf: Histogram of raw p-values.
* hist_p_values_cutoff_Pos.pdf: Histogram of p-values less than 0.05.
* heatmap_CNA_Pos_df2.pdf: Heatmap of CNA-TCs.
* 3heatmap_TACNA_CNA_TCs.pdf: Heatmap of TACNA-CNA-TCs.
* immune_cna_tc_for_projection_spatialsamples.tsv: Merged dataframe of unique TCs.


#######################################################################################

Script: merge_immune_tcs_gobp_reactome.R

Description:
This script processes immune-TCs identified from Reactome and GOBP databases. It combines and processes these TCs, and generates heatmaps for visual analysis. The primary goal is to identify and visualize significant TCs related to immune response.

Data Processing Steps:

* Extract Column Names: Extract column names from GO-BP and Reactome gene set files.
* Create DataFrame of Unique TCs: Combine unique TCs, check their presence in GO-BP and Reactome datasets, and save the dataframe.
* Subset ICA Flipped Data: Load ICA results, subset to keep only immune TCs, and save the subsetted dataframe.
* Process GO-BP Data: Subset GO-BP gene sets based on overlaps, process gene set names, further subset to match gene sets, and save the dataframe.
* Process Reactome Data: Subset Reactome gene sets based on overlaps, process gene set names, further subset to match gene sets, and save the dataframe.
* Generate Heatmaps: Define color palette, shorten names for readability, generate and save heatmaps for both Reactome and GO-BP data, and center the scale to zero for additional heatmaps.

Input Files:
* Overlapping gene sets from GO-BP (tab-separated file).
* Overlapping gene sets from Reactome (tab-separated file).
* ICA (Independent Component Analysis) results (tab-separated file).
* GOBP gene sets with Bonferroni corrected p-values (dataframe).
* Reactome gene sets with Bonferroni corrected p-values (dataframe).

Output Files:
* Presence of TCs in GOBP and Reactome datasets (tab-separated file).
* Subset of ICA results corresponding to immune TCs (tab-separated file).
* Subset of GOBP and Reactome gene sets corresponding to immune TCs (tab-separated files).
* Heatmaps for immune Reactome and GOBP data (PDF files).


#######################################################################################

Script: genomic_mapping.R

Description:
This script generates genomic ordered plots of CNA-TCs derived from the TACNA method. The plots visualize the genomic distribution of TC weights across different chromosomes, highlighting their genomic mapping.

Data Processing Steps:
* Set Up: Set the working directory and define parameters.
* Load Data: Read gene annotation, CNA ICA results, and columns to plot files.
* Preprocess Data:
o Filter and compute middle base pair positions for genes.
o Merge gene annotations with CNA results.
o Convert chromosome annotations to numeric values.
* Prepare Plotting Data: Order data by chromosome and base pair, and create chromosome annotations.
* Generate Plots: Define and use a function to create and save genomic ordered plots for each TC.
This script provides a visual representation of the genomic distribution of TCs derived from CNA analysis, aiding in the identification of significant genomic regions associated with these components.

Input Files:
* Gene Annotation File: Mapping of genes to their chromosomal locations (tab-separated file).
* CNA ICA Results File: ICA results for CNAs (tab-separated file).
* Columns to Plot File: List of TCs to be plotted (tab-separated file).

Output Files:
Genomic Ordered Plot Files: PNG files of genomic ordered plots for each TC.


#######################################################################################

Script: plotting_tcs_on_spatial_samples.R

Description:
This script is designed to plot TCs on spatial transcriptomics samples. The script reads in mixing matrices from multiple datasets, processes the data, and generates spatial plots for the specified TCs. The generated plots are then saved as a PDF file.

Steps of Data Processing:
* Extract Column Names
* Read and combine Mixing Matrices
* Clean Workspace
* Read Spatial Data
* Generate Plots
    Define functions to plot raw images and TCs.
    Generate and save plots for each sample and TC combination.

Input Files:
* Mixing Matrices
TSV (Tab-Separated Values): These files contain the mixing matrices from different datasets (BC1 to BC5). Each file has rows representing samples and columns representing independent components.
* Spatial Data
RDS (R Data File): A single RDS file containing spatial transcriptomics data for all samples. This file includes the spatial coordinates and associated metadata.

Output Files:
PDF file containing spatial plots of TCs for each sample.


#######################################################################################

Script: plotting_CNA_TCs_on_spatial_samples.R

Description:
This script visualizes TCs in spatial transcriptomics samples. It identifies CNA-TCs with 50 or more genes, loads mixing matrices from multiple spatial samples, filters the TCs, and plots them on spatial sample images. The visualizations help understand the spatial distribution of TCs in various cancer samples.

Data Processing Steps:
* Load Data: Read CNA-TCs data and identify unique TCs.
* Filter TCs: Select TCs with 50 or more genes.
* Load Mixing Matrices: Read and combine mixing matrices from multiple samples.
* Filter Mixing Matrices: Filter to include only the selected TCs.
* Load Spatial Data: Read merged spatial data.
* Generate Plots: Plot the selected TCs on spatial sample images.
* Save Plots: Save the plots as a PDF file.

Input Files:
* CNA-TCs Data: Tab-separated values file containing CNA-TC data.
* Spatial Data: RDS file containing merged spatial data for all Visium samples.
* Mixing Matrices: Multiple tab-separated files containing mixing matrices for various samples.

Output Files:
* Filtered TCs: Tab-separated values file of TCs with 50 or more genes.
* Combined Mixing Matrix: Tab-separated values file of the combined and filtered mixing matrix.
* TC Plots: PDF file containing plotted TCs for each spatial sample.


####################################################################################### 

Script: check_colocalization_ks_slurm_part1.R

Description
This script performs colocalization correlation analysis between TCs in spatial transcriptomics samples. It uses kernel density estimates (KDE) to visualize spatial distributions, computes colocalization correlations, and calculates coloc scores to measure the interaction between TCs.

Steps:
* Load Libraries and Set Up: Load necessary libraries and set the number of cores for parallel processing.
* Read Input Files: Load the corrected mixing matrix and merged spatial data.
* Create Mapping: Generate a mapping matrix for each image, aligning spatial coordinates with the mixing matrix.
* KDE Calculation: Calculate KDE estimates for each TC across images, using a quantile cutoff to create density masks.
* Save KDE Estimates: Save the KDE estimates as RDS files.
* Colocalization Correlation Calculation: Compute Pearson correlation coefficients between KDE estimates of TCs.
* Save Correlation Results: Save the correlation results in a tab-separated values file.
* Calculate Coloc Scores: Compute coloc scores from the correlation data, summarizing the interaction between TCs.
* Save Coloc Scores: Save the coloc scores in a tab-separated values file.

Input Files:
* Corrected Mixing Matrix: corrected_mix_mat.tsv
* Merged Visium Spatial Data: ALL_VISIUM_FILES_MERGED.rds

Output Files:
* KDE Estimates: RDS files containing KDE estimates for TCs.
* Colocalization Correlations: Tab-separated values file of the Pearson correlation coefficients between TCs.
* Coloc Scores: Tab-separated values file of the calculated coloc scores.


#######################################################################################

Script: check_colocalization_ks_slurm_part2.R

Description:
This script performs colocalization analysis between TCs derived from independent component analysis (ICA) on spatial transcriptomics data. It calculates KDE and performs correlation analysis to identify patterns of colocalization between different TCs in spatial data from Visium slides.

Steps of Data Processing:
* Setup and Library Loading: Load necessary R libraries and set up parallel processing.
* Read Input Data: Load the corrected mixture matrix and Visium spatial transcriptomics data.
* Define Helper Functions: Define functions for creating mappings between spatial coordinates and TCs, and for applying quantile cutoffs.
* Calculate KDE Estimates: For each image and TC, calculate KDE estimates and save the results as PNG plots and RDS files.
* Calculate Colocalization Correlation: Compute Pearson correlation coefficients between KDE estimates of TC pairs, and save the results as a TSV file.
* Calculate Colocalization Scores: Compute colocalization scores based on the correlation results and save individual and combined TSV files for further analysis.

Input Files:
* Corrected mixture matrix (TSV file) containing samples as rows and TCs as columns.
* Merged Visium spatial transcriptomics data (RDS file).

Output Files:
* KDE Estimates:
* Directory containing PNG images of KDE plots for each TC and image.
* RDS file containing KDE estimates for TCs.
* Colocalization Correlation:
* TSV file containing Pearson correlation coefficients and p-values for all TC pairs.
* Colocalization Score:
* TSV files for each image, containing colocalization scores.
* Combined TSV file containing colocalization scores for all images.



#######################################################################################

Script: check_colocalization_ks_slurm_part3.R

Description:
This script performs colocalization analysis between TCs derived from independent component analysis (ICA) on spatial transcriptomics data. It calculates KDE and performs correlation analysis to identify patterns of colocalization between different TCs in spatial data from Visium slides.

Steps of Data Processing:
Setup
* Load libraries and custom functions.
* Set up parallel processing and output directory.
Read Data
* Load expression matrix and Visium data.
Kernel Density Estimation (KDE)
* Calculate KDE for each cell type in each image.
* Save and plot KDE results.
Colocalization Correlation
* Calculate Pearson correlations for cell type combinations.
* Save correlation results.
Colocalization Scores
* Filter significant correlations.
* Compute and save colocalization scores for each image.

Input Files:
* Corrected mixture matrix (TSV file) containing samples as rows and TCs as columns.
* Merged Visium spatial transcriptomics data (RDS file).

Output Files:
* KDE Estimates:
* Directory containing PNG images of KDE plots for each TC and image.
* RDS file containing KDE estimates for TCs.
* Colocalization Correlation:
* TSV file containing Pearson correlation coefficients and p-values for all TC pairs.
* Colocalization Score:
* TSV files for each image, containing colocalization scores.
* Combined TSV file containing colocalization scores for all images.



#######################################################################################

Script: identify_correlated_TCs.R

Description:
This script identifies TCs after colocalization analysis by evaluating the correlation values between pairs of TCs on spatial samples. It considers three types of TC pairs: CNA-TC + CNA-TC, CNA-TC + Immune-TC, and Immune-TC + Immune-TC. Heatmaps are generated for each sample showing the correlated TCs.

Data Processing Steps:
* Load Data: Read in colocalization scores, immune TCs, and CNA TCs.
* Identify Common TCs: Determine common TCs between immune and CNA sets.
* Rename Columns: Standardize TC column names.
* Subset Data: Identify pairs of TCs present in colocalization scores for both immune and CNA sets.
* Generate Heatmaps: Create heatmaps for each sample based on colocalization scores.
* Identify Top TC Pairs: Find the highest and lowest correlation values among CNA TCs, immune TCs, and mixed CNA + immune TCs.

Input Files:
* Colocalization Scores: spatial_kernels_coloc_score4_cutoff2_all.tsv
* Immune TCs: immune_ica_flipped_df.tsv
* CNA TCs: cna_ica_flipped_df.tsv
* All TCs Combined: immune_cna_tc_for_projection_spatialsamples.tsv

Output Files:
* Matching CNA Pairs: matching_cna_pairs.tsv
* Matching Immune Pairs: matching_immune_pairs.tsv
* CNA Heatmap: 4colocalization_cna_heatmap.pdf
* Immune Heatmap: 4colocalization_immune_heatmap.pdf
* Unique TC Colocalization (CNA): unique_tc_coloc_cna.tsv
* Unique TC Colocalization (Immune): unique_tc_coloc_immune.tsv
* Subset Colocalization Scores: cnatc1_immunetc2_subset_coloc_scores.tsv



#######################################################################################

Script: plotting_top_colocalized_tcs_on_spatial_samples.R

Description:
This script helps in visualizing the spatial distribution of TCs, which is crucial for understanding their colocalization and correlation in different spatial transcriptomics samples, focusing on those with the highest and lowest correlation values.

Data Processing Steps:
* Load necessary libraries and set the output directory.
* Read the correlation files for CNA-TCs and immune-TCs.
* Define the highest and lowest correlation values for different TC pairs.
* Load spatial data and results from the projection of all TCs on all spots.
* Define TCs to plot.
* Define plot functions to visualize raw and specific TC images on spatial data.
* Iterate through all specified images and TCs, generating plots for each.
* Save the combined plots to a PDF file.

Input Files:
Correlation data for CNA-TCs and immune-TCs.
Spatial transcriptomics data.
Results from the projection of TCs on all spots in all samples. 

Output Files:
PDF file containing plots of the selected TCs on the specified spatial samples.


#######################################################################################

Script: TC_correlation_histograms.R

Description: 
This R script generates histograms and bar plots of colocalization correlation values for TCs across spatial transcriptomics cancer samples. It visualizes maximum correlation scores for CNA-TCs, immune-TCs, and immune-CNA-TCs collectively and individually. These visualizations help analyze the distribution and relationships of TC colocalization scores across cancer samples.

Data Processing Steps:
* Load correlation data for CNA-TCs, immune-TCs, and immune-CNA-TCs from respective input files.
* Calculate maximum correlation scores for each cancer sample.
* Generate bar plots for maximum correlation scores and histograms for correlation distributions across cancer samples using base R functions and reshape2 library for data manipulation. Adjust plot parameters for layout and visualization clarity.
* Save the generated plots as PDF files in the specified output directory.

Input Files:
* unique_tc_coloc_cna.tsv: Tab-separated file containing correlation scores of CNA-TCs.
* unique_tc_coloc_immune.tsv: Tab-separated file containing correlation scores of immune-TCs.
* cnatc1_immunetc2_subset_coloc_scores.tsv: Tab-separated file containing subset correlation scores of immune-CNA-TCs.

Output Files:
* Pdf file of bar plot showing maximum correlation scores of CNA-TCs across cancer samples.
* Pdf file of histograms showing correlation distributions of CNA-TCs for each cancer sample.
* Pdf file of bar plot showing maximum correlation scores of immune-TCs across cancer samples.
* Pdf file of histograms showing correlation distributions of immune-TCs for each cancer sample.
* Pdf file of bar plot showing maximum correlation scores of immune-CNA-TCs across cancer samples.
* Pdf file of histograms showing correlation distributions of immune-CNA-TCs for each cancer sample.


#######################################################################################

Script: subset_old_corrected_mixing_matrix.R

Description:
This script subsets the previous-wrong corrected mixing matrix (derived from wrong analysis) by retaining only the 476 CNA-TCs. It aligns and saves the subsetted matrix based on the intersection with a previously defined set of CNA-TCs, facilitating downstream analysis.

Input files:
Tab-separated files (.tsv): contains genomic data of CNA-TCs, immune-TCs, and all TCs for spatial sample projection, as well as an old-wrong corrected mixing matrix.

Output files:
Tab-separated file (.tsv): Subsetted corrected mixing matrix containing only CNA-TCs, optimized for subsequent data processing.


#######################################################################################

Script: creating_new_mixing_matrix.R

Description:
This R script creates a new corrected mixing matrix by combining CNA-TCs and the new, correct immune-TCs from various cancer spatial transcriptomics samples. It loads these data files, merges multiple mixing matrices, renames columns appropriately, removes duplicates, and saves the final matrix, which is the new corrected mixing matrix with the old, CNA TCs. 

Steps of Data Processing:
* Load Data: Load CNA-TCs, immune-TCs, and various mixing matrices from specified directories.
* Combine Mixing Matrices: Combine all mixing matrices into a single dataframe.
* Rename Columns: Rename columns in the combined mixing matrix to follow a standard pattern.
* Merge Matrices: Merge the corrected CNA-TCs mixing matrix with the combined mixing matrix.
* Remove Duplicates: Identify and remove duplicate column names.
* Save and Read: Save the final corrected mixing matrix and read it back for further analysis.
* Spatial Data: Read in spatial transcriptomics data from a specified RDS file for additional processing.

Input Files:
* Tab-delimited file (.tsv): Contains corrected mixing matrix for CNA-TCs.
* Tab-delimited file (.tsv): Contains all TCs data for spatial transcriptomics analysis.
* Tab-delimited file (.tsv): Contains correct immune-TCs data.

Output Files:
Tab-delimited file (.tsv): Contains the new corrected mixing matrix.


#######################################################################################

Script: cna_tcs_heatmap_chr_annotations.R

Description:
This R script is to visualize correlated CAN-TCs data across spatial transcriptomics samples, generating heatmaps with chromosome-scale annotations.

Data Processing Steps:
* Load Libraries: Load necessary R libraries (Seurat, dplyr, pheatmap, etc.).
* Set Working Directory: Set directory to where input files are located.
* Load Data: Read input files with matched CNA-TC pairs and chromosome mappings.
* Add Chromosome Annotations: Merge chromosome mappings into the CNA-TC pairs data.
* Generate Heatmaps:
* Loop through unique identifiers.
* Subset data, reshape, and prepare chromosome annotations.
* Plot heatmaps with clustering and chromosome annotations using pheatmap.
* Save resulting heatmaps as PDF files.

Input Files:
* Tab-separated file: Contains matched pairs of CNAs and TCs with correlation scores.
* Tab-separated file: Includes chromosome mappings and TC identifiers.

Output Files:
* PDF file: Heatmaps of correlated CNA-TCs using rainbow chromosome annotation.
* PDF file: Heatmaps using distinct Polychrome color palette for chromosome annotation.
* PDF file: Corrected chromosome order and rainbow chromosome annotation.


#######################################################################################

Script:  explore_top_colocalized_TCs

Description:
This R script is for exploring the biological information of the top colocalized CNA- and immune-TCs. 

Data Processing Steps:
* Load Libraries: Import necessary R libraries (dplyr, data.table).
* Load Data: Read various input files containing immune-TCs, CNA-TCs, positional data, enrichment matrices, ICA results, EntrezID mappings, and colocalization scores.
* Explore Top Pathways and TCs:
* Identify top pathways associated with immune-TCs using GOBP and Reactome databases.
* Investigate gene weights and corresponding EntrezIDs for top TCs.
* Examine colocalization scores for specific TC pairs across cancer samples.
* Gene and Pathway Analysis:
o Retrieve gene names corresponding to EntrezIDs for top TCs.
o Analyze correlation scores and pathway enrichments using heatmap examples.

Input Files:
* Tab-separated file: Contains immune TC data for pathway enrichment and analysis.
* Tab-separated file: Includes CNA-TC data for correlation and pathway exploration.
* Tab-separated file: Merged dataset of unique immune and CNA TCs.
* Tab-separated file: Positional data for cytoband exploration of CNA-TCs.
* Tab-separated file: Enrichment matrix from GOBP database.
* Tab-separated file: Enrichment matrix from Reactome database.
* Tab-separated file: ICA results for top gene weights per pathway.
* Tab-separated file: EntrezID mapping for gene name lookup.
* Tab-separated file: Colocalization scores of CNA-TCs and immune-TCs across cancer samples.

Output Files:
* Analysis results and explorations for top immune and CNA TCs.
* Pathway annotations and gene name mappings.
* Colocalization scores for specific TC pairs across cancer samples.









References: 
OpenAI. (2024). Script Summary and Optimization by ChatGPT (Version 2.5): Assisted in script summary of the R script provided by user, with added comments for clarity and optimizations. 



