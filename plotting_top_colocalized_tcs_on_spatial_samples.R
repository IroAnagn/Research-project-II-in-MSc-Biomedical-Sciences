############## Plotting the top colocalized TCs on spatial transcriptomics samples ####################



# Source additional functions from GitHub
source("https://raw.githubusercontent.com/loipf/coding_snippets/master/R/small_functions.R")

# Load necessary libraries
library(RColorBrewer)
library(Seurat)
library(circlize)
library(data.table)
library(ggplot2)
library(cowplot)

# Set output directory
OUTPUT_DIR = '/home/iro/Colocalization_results/correlation plots'
setwd(OUTPUT_DIR)

################ Load the correlation files ################
unique_tc_coloc_cna <- read.table("/home/iro/Colocalization_results/unique_tc_coloc_cna.tsv", sep = '\t', header = TRUE)
unique_tc_coloc_immune <- read.table("/home/iro/Colocalization_results/unique_tc_coloc_immune.tsv", sep = '\t', header = TRUE)
tc_coloc_cna_immune <- read.table("/home/iro/Colocalization_results/cnatc1_immunetc2_subset_coloc_scores.tsv", sep = '\t', header = TRUE)

## Highest and lowest correlation values
## Highest correlation value between CNA-TCs: TC3557_TC6547 in OC2.
## Lowest correlation value between CNA-TCs: TC8914_TC9155 in CRC2. 
## Highest correlation value between immune-TCs: TC372_TC641  in BC3. 
## Lowest correlation value between immune-TCs: TC4121_TC641  in GB1.  
## Highest correlation value for CNA TC + Immune TC: TC8628_TC8538  in BC1.   
## Lowest correlation value for CNA TC + Immune TC: TC7834_TC8928  in OC2. 
## No correlation - zero value: TC3567_TC6018 in BC1.  
## No correlation - zero value: TC9517_TC4669 in BC5. 

################ Plot these TCs on all spatial samples ##########
# Read the spatial data
visium <- readRDS("/home/arkajyotibhattacharya/Projects/spatial_transcriptomics_marco/Data/ALL_VISIUM_FILES_MERGED.rds")

# Read the file containing the results from projection of all TCs on all spots in all samples
corrected_mix_mat <- read.table("/home/iro/Colocalization_results/new_corrected_mixing_matrix.tsv", sep = '\t', header = TRUE)

# Select TCs to plot
tcs_to_plot <- c("TC4121", "TC641")

### Plot functions
plot_raw_image <- function(image_id, corrected_mix_mat, visium) {
  mapping_file <- visium@images[[image_id]]@coordinates
  
  # Adjust mapping for specific image if needed
  if(image_id == "OC1") {
    mapping_file <- subset(mapping_file, col < 185)
  }
  
  common_samples <- intersect(rownames(mapping_file), rownames(corrected_mix_mat))
  mapping_file <- mapping_file[common_samples,]
  corrected_mix_mat_current <- corrected_mix_mat[common_samples,]
  corrected_mix_mat_transposed <- t(corrected_mix_mat_current)
  
  new.seurat.object <- CreateSeuratObject(counts = corrected_mix_mat_transposed, assay = "Spatial")
  new.seurat.object@images$image <- new(
    Class = 'VisiumV1',
    assay = "Spatial",
    key = "image_",
    coordinates = mapping_file,
    image = visium@images[[image_id]]@image,
    scale.factors = visium@images[[image_id]]@scale.factors,
    spot.radius = visium@images[[image_id]]@spot.radius
  )
  
  TCs_to_plot_index <- colnames(corrected_mix_mat)[1]
  p <- SpatialFeaturePlot(new.seurat.object, features = TCs_to_plot_index, slot = "counts", stroke = NA, image.alpha = 1, pt.size.factor = 1.8, alpha = 0, interactive = FALSE) & NoLegend() & NoAxes()
  p
}

plot_tc_image <- function(image_id, tc_to_plot, corrected_mix_mat, visium) {
  mapping_file <- visium@images[[image_id]]@coordinates
  
  if(image_id == "OC1") {
    mapping_file <- subset(mapping_file, col < 185)
  }
  
  common_samples <- intersect(rownames(mapping_file), rownames(corrected_mix_mat))
  mapping_file <- mapping_file[common_samples,]
  corrected_mix_mat_current <- corrected_mix_mat[common_samples,]
  corrected_mix_mat_transposed <- t(corrected_mix_mat_current)
  
  scale_fill_custom <- ggplot2::scale_fill_gradientn(limits = c(-6, 6), breaks = seq(-6, 6, by = 3), colours = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100), oob = scales::squish)
  
  new.seurat.object <- CreateSeuratObject(counts = corrected_mix_mat_transposed, assay = "Spatial")
  new.seurat.object@images$image <- new(
    Class = 'VisiumV1',
    assay = "spatial",
    key = "image_",
    coordinates = mapping_file,
    image = visium@images[[image_id]]@image,
    scale.factors = visium@images[[image_id]]@scale.factors,
    spot.radius = visium@images[[image_id]]@spot.radius
  )
  
  p <- SpatialFeaturePlot(new.seurat.object, features = tc_to_plot, slot = "counts", stroke = NA, image.alpha = 0.2, alpha = 1, interactive = FALSE) & scale_fill_custom & NoAxes() & NoLegend() & ggtitle(tc_to_plot)
  p
}

### Iterate through all images
images_to_plot <- c("GB1")

plot_list_tcs <- lapply(images_to_plot, function(curr_img) {
  print(paste0("### start ", curr_img))
  plot_list_tcs_curr_tc <- lapply(c("text", 'raw', tcs_to_plot), function(curr_tc) {
    print(paste0("### ", curr_tc))
    p <- NULL
    if(curr_tc == 'text') {
      p <- ggplot() + geom_text(aes(label = I(curr_img), y = 0.5, x = 0.2), hjust = 0) + xlim(0, 1) + ylim(0, 1) + theme_void()
    } else if(curr_tc == 'raw') {
      p <- plot_raw_image(curr_img, corrected_mix_mat, visium) + theme(plot.margin = margin(l = -0.5, unit = "in"))
    } else {
      p <- plot_tc_image(curr_img, curr_tc, corrected_mix_mat, visium) + theme(plot.margin = margin(l = -0.5, unit = "in"))
    }
    p
  })
  plot_list_tcs_curr_tc
})

plot_list_tcs <- unlist(plot_list_tcs, recursive = FALSE)

pdf('TC4121_TC641_in_GB1.pdf', width = 3 * length(tcs_to_plot), height = 2 * length(images_to_plot))
plot_grid(plotlist = plot_list_tcs, nrow = length(images_to_plot), rel_widths = c(0.6, rep(1, length(tcs_to_plot) + 2)), scale = 1)
dev.off()
