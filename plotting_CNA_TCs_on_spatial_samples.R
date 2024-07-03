############## Plotting TCs on spatial samples ####################
################ Transformation to S4 objects #########################

source("https://raw.githubusercontent.com/loipf/coding_snippets/master/R/small_functions.R")

# Load necessary libraries
library(RColorBrewer)
library(Seurat)
library(circlize)
library(data.table)
library(ggplot2)
library(cowplot)

# Set output directory
OUTPUT_DIR = '/home/iro/Projection_results_wrong/Projection plots'
setwd(OUTPUT_DIR)


################ load the CNA-TCs file  ################
# Read the CNA-TCs file into a dataframe
CNA_TCs <- data.frame(read.table("/home/iro/GSEA_results/TACNA_CNA_TCs.tsv", sep = '\t', header = TRUE))

unique_TCs <- table(CNA_TCs$Component)
unique_TCs_df <- as.data.frame(unique_TCs)

############ Check which TCs have 50 or more genes ############
cna_tcs_50genes <- unique_TCs_df$Freq >= 50
sum(cna_tcs_50genes)

CNA_50genes_df <- unique_TCs_df[cna_tcs_50genes, ]

# save it 
write.table(CNA_50genes_df, file = "/home/iro/Projection_results_wrong/CNA_TCs_50genes.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
# Continue from this step: 
CNA_50genes_df <- read.table(file = "/home/iro/Projection_results_wrong/CNA_TCs_50genes.tsv", sep = '\t', header = TRUE)


########### plot on all spatial samples only these 25 CNA-TCs ##########

############# read in all mixing matrices from projection 
mix_mat_bc1 = data.frame(fread("/home/iro/Projection_results/BC1/Project_data_on_independent_components_{865ba4ef-97fe-4299-bd00-c0f4f7fa2571}/mixing_matrix.tsv"), row.names = 1)
mix_mat_bc2 = data.frame(fread("//home/iro/Projection_results/BC2/Project_data_on_independent_components_{ac65f7d4-6fbc-4de3-b046-50eab769705f}/mixing_matrix.tsv"), row.names = 1)
mix_mat_bc3 = data.frame(fread("/home/iro/Projection_results/BC3/Project_data_on_independent_components_{eade2a1a-9903-49fc-864b-bf0dceca5c17}/mixing_matrix.tsv"), row.names = 1)
mix_mat_bc4 = data.frame(fread("/home/iro/Projection_results/BC4/Project_data_on_independent_components_{63ab3898-a12d-41c6-a4ba-679f91bf6884}/mixing_matrix.tsv"), row.names = 1)
mix_mat_bc5 = data.frame(fread("/home/iro/Projection_results/BC5/Project_data_on_independent_components_{c83598c9-4fd3-4f99-980a-4cf6c34ba7fb}/mixing_matrix.tsv"), row.names = 1)
mix_mat_cc1 = data.frame(fread("/home/iro/Projection_results/CC1/Project_data_on_independent_components_{7b84a2d5-cc65-4eca-97c5-a3bc5db4adeb}/mixing_matrix.tsv"), row.names = 1)
mix_mat_crc1 = data.frame(fread("/home/iro/Projection_results/CRC1/Project_data_on_independent_components_{67c0032f-e417-4bef-ae5c-62a264f3537f}/mixing_matrix.tsv"), row.names = 1)
mix_mat_crc2 = data.frame(fread("/home/iro/Projection_results/CRC2/Project_data_on_independent_components_{f333622c-9eac-42e1-908f-40982776efbc}/mixing_matrix.tsv"), row.names = 1)
mix_mat_crc3 = data.frame(fread("/home/iro/Projection_results/CRC3/Project_data_on_independent_components_{d19e68df-738d-4296-ae8e-540a5b54191a}/mixing_matrix.tsv"), row.names = 1)
mix_mat_gb1 = data.frame(fread("/home/iro/Projection_results/GB1/Project_data_on_independent_components_{b5fa7284-c066-4467-9a23-e0fa81953aa1}/mixing_matrix.tsv"), row.names = 1)
mix_mat_gb2 = data.frame(fread("/home/iro/Projection_results/GB2/Project_data_on_independent_components_{321125e6-5db0-4e77-8297-8d2042020608}/mixing_matrix.tsv"), row.names = 1)
mix_mat_idc1 = data.frame(fread("/home/iro/Projection_results/IDC1/Project_data_on_independent_components_{21086816-b553-4f71-b7e5-df4de2903d29}/mixing_matrix.tsv"), row.names = 1)
mix_mat_lc1 = data.frame(fread("/home/iro/Projection_results/LC1/Project_data_on_independent_components_{9088d440-4ce3-4a8d-8d70-67e4a09f0413}/mixing_matrix.tsv"), row.names = 1)
mix_mat_lc2 = data.frame(fread("/home/iro/Projection_results/LC2/Project_data_on_independent_components_{68d8a498-04ef-487a-bead-03804ba6c400}/mixing_matrix.tsv"), row.names = 1)
mix_mat_m1 = data.frame(fread("/home/iro/Projection_results/M1/Project_data_on_independent_components_{633c35bd-333d-42cb-a4dd-3271eeb9961a}/mixing_matrix.tsv"), row.names = 1)
mix_mat_oc1 = data.frame(fread("/home/iro/Projection_results/OC1/Project_data_on_independent_components_{0d99cfb4-e7f5-42bf-b865-86fc3cc8e140}/mixing_matrix.tsv"), row.names = 1)
mix_mat_oc2 = data.frame(fread("/home/iro/Projection_results/OC2/Project_data_on_independent_components_{52e515d9-9cac-42ed-8701-cf292ec1b47b}/mixing_matrix.tsv"), row.names = 1)
mix_mat_oc3 = data.frame(fread("/home/iro/Projection_results/OC3/Project_data_on_independent_components_{4ef547f1-2aa7-485f-ab44-2ac483174774}/mixing_matrix.tsv"), row.names = 1)
mix_mat_pc1 = data.frame(fread("/home/iro/Projection_results/PC1/Project_data_on_independent_components_{85773807-e91c-4e21-80d9-fe82b35aee23}/mixing_matrix.tsv"), row.names = 1)
mix_mat_pc2 = data.frame(fread("/home/iro/Projection_results/PC2/Project_data_on_independent_components_{a26d872a-61d1-41d4-9ee2-e50c2c4a71c8}/mixing_matrix.tsv"), row.names = 1)
mix_mat_pc3 = data.frame(fread("/home/iro/Projection_results/PC3/Project_data_on_independent_components_{e41911f5-f0a2-420f-a8e3-33452dab9ada}/mixing_matrix.tsv"), row.names = 1)

  
# Combine the mixing matrices into one dataframe
corrected_mix_mat = as.data.frame(rbind(mix_mat_bc1, mix_mat_bc2))
corrected_mix_mat = as.data.frame(rbind(corrected_mix_mat, mix_mat_bc3))
corrected_mix_mat = as.data.frame(rbind(corrected_mix_mat, mix_mat_bc4))
corrected_mix_mat = as.data.frame(rbind(corrected_mix_mat, mix_mat_bc5))
corrected_mix_mat = as.data.frame(rbind(corrected_mix_mat, mix_mat_cc1))
corrected_mix_mat = as.data.frame(rbind(corrected_mix_mat, mix_mat_crc1))
corrected_mix_mat = as.data.frame(rbind(corrected_mix_mat, mix_mat_crc2))
corrected_mix_mat = as.data.frame(rbind(corrected_mix_mat, mix_mat_crc3))
corrected_mix_mat = as.data.frame(rbind(corrected_mix_mat, mix_mat_gb1))
corrected_mix_mat = as.data.frame(rbind(corrected_mix_mat, mix_mat_gb2))
corrected_mix_mat = as.data.frame(rbind(corrected_mix_mat, mix_mat_idc1))
corrected_mix_mat = as.data.frame(rbind(corrected_mix_mat, mix_mat_lc1))
corrected_mix_mat = as.data.frame(rbind(corrected_mix_mat, mix_mat_lc2))
corrected_mix_mat = as.data.frame(rbind(corrected_mix_mat, mix_mat_m1))
corrected_mix_mat = as.data.frame(rbind(corrected_mix_mat, mix_mat_oc1))
corrected_mix_mat = as.data.frame(rbind(corrected_mix_mat, mix_mat_oc2))
corrected_mix_mat = as.data.frame(rbind(corrected_mix_mat, mix_mat_oc3))
corrected_mix_mat = as.data.frame(rbind(corrected_mix_mat, mix_mat_pc1))
corrected_mix_mat = as.data.frame(rbind(corrected_mix_mat, mix_mat_pc2))
corrected_mix_mat = as.data.frame(rbind(corrected_mix_mat, mix_mat_pc3))


# Remove redundant objects from the workspace
rm(mix_mat_bc2)
rm(mix_mat_bc3)
rm(mix_mat_bc4)
rm(mix_mat_bc5)
rm(mix_mat_cc1)
rm(mix_mat_crc1)
rm(mix_mat_crc2)
rm(mix_mat_crc3)
rm(mix_mat_gb1)
rm(mix_mat_gb2)
rm(mix_mat_idc1)
rm(mix_mat_lc1)
rm(mix_mat_lc2)
rm(mix_mat_m1)
rm(mix_mat_oc1)
rm(mix_mat_oc2)
rm(mix_mat_oc3)
rm(mix_mat_pc1)
rm(mix_mat_pc2)
rm(mix_mat_pc3)


#	Read the spatial data 
visium = readRDS("/home/arkajyotibhattacharya/Projects/spatial_transcriptomics_marco/Data/ALL_VISIUM_FILES_MERGED.rds")

colnames(corrected_mix_mat) = gsub(pattern = "consensus.independent.component.", replacement = "TC", x = colnames(corrected_mix_mat))

### find the 25 CNA-TCs in corrected_mix_mat when compared with CNA_50genes_df
# Extract TC numbers from CNA_50genes_df$Var1
tc_numbers <- as.character(CNA_50genes_df$Var1)

# Convert TC numbers to the format used in corrected_mix_mat
tc_names <- paste0("TC", tc_numbers)

# Check which TCs exist in corrected_mix_mat
existing_tcs <- tc_names[tc_names %in% colnames(corrected_mix_mat)]

write.table(corrected_mix_mat, file = "/home/iro/Projection_results_wrong/corrected_mix_mat.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)

corrected_mix_mat <- read.table(file = "/home/iro/Projection_results_wrong/corrected_mix_mat.tsv", sep = '\t', header = TRUE)


# Filter corrected_mix_mat to include only the 25 TCs
CNA_corrected_mix_mat <- corrected_mix_mat[, existing_tcs, drop = FALSE]


# save it 
write.table(CNA_corrected_mix_mat, file = "/home/iro/Projection_results_wrong/CNA_corrected_mix_mat.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
# Continue from this step: 
CNA_corrected_mix_mat <- read.table(file = "/home/iro/Projection_results_wrong/CNA_corrected_mix_mat.tsv", sep = '\t', header = TRUE)



###### select here which TCs to plot
# tcs_to_plot = colnames(corrected_mix_mat)[1:20]
tcs_to_plot = colnames(CNA_corrected_mix_mat)


############################################

### plot functions


plot_raw_image = function(image_id, CNA_corrected_mix_mat = CNA_corrected_mix_mat, visium = visium) {
  
  mapping_file <- visium@images[[image_id]]@coordinates
  
  if(image_id =="OC1") {
    
    mapping_file = subset(mapping_file, col<185)  ### outside image region
    
  }
  
  
  common_samples = intersect(rownames(mapping_file),rownames(CNA_corrected_mix_mat) )
  
  mapping_file = mapping_file[common_samples,]
  
  corrected_mix_mat_current = CNA_corrected_mix_mat[common_samples,]
  
  corrected_mix_mat_transposed = t(corrected_mix_mat_current)
  
  
  
  new.seurat.object = CreateSeuratObject(counts = corrected_mix_mat_transposed, assay = "Spatial" )
  
  new.seurat.object@images$image = new(
    
    Class = 'VisiumV1'
    
    ,assay = "Spatial"
    
    ,key = "image_"
    
    ,coordinates = mapping_file
    
    ,image = visium@images[[image_id]]@image
    
    ,scale.factors = visium@images[[image_id]]@scale.factors
    
    ,spot.radius = visium@images[[image_id]]@spot.radius
    
  )
  
  
  TCs_to_plot_index = colnames(CNA_corrected_mix_mat)[1]
  
  p = SpatialFeaturePlot(new.seurat.object
                         
                         , features =TCs_to_plot_index
                         , slot = "counts"
                         
                         , stroke = NA
                         
                         , image.alpha = 1
                         
                         , pt.size.factor = 1.8
                         
                         , alpha = 0
                         
                         # , interactive = FALSE)&NoLegend()&NoAxes()
                         
                         , interactive = FALSE)&NoAxes()
  
  p
  
}



plot_tc_image = function(image_id, tc_to_plot, CNA_corrected_mix_mat = CNA_corrected_mix_mat, visium = visium) {
  
  mapping_file <- visium@images[[image_id]]@coordinates
  
  
  
  if(image_id =="OC1") {
    
    mapping_file = subset(mapping_file, col<185)  ### outside image region
    
  }
  
  
  
  common_samples = intersect(rownames(mapping_file),rownames(CNA_corrected_mix_mat) )
  
  mapping_file = mapping_file[common_samples,]
  
  corrected_mix_mat_current = CNA_corrected_mix_mat[common_samples,]
  
  
  
  # corrected_mix_mat_current$ = range01(corrected_mix_mat_current$) ## cna per image
  
  corrected_mix_mat_transposed = t(corrected_mix_mat_current)
  
  
  if(tc_to_plot==""){
    
    scale_fill_custom = ggplot2::scale_fill_gradient(limits = c(0,1), low = "white", high = "red" )
    
  } else {
    
    # scale_fill_custom = ggplot2::scale_fill_gradient2(limits = c(-20,20), breaks=seq(-20, 20, by = 2), low = "#00006b", mid = "white", high = "#a70000", midpoint = 0, oob = scales::squish)
    
    # scale_fill_custom = ggplot2::scale_fill_gradient2(limits = c(-10,10), breaks=seq(-10, 10, by = 2), low = "#00006b", mid = "white", high = "#a70000", midpoint = 0, oob = scales::squish)
    
    # scale_fill_custom = ggplot2::scale_fill_gradient2(limits = c(-15,15), breaks=seq(-15, 15, by = 2), low = "#00006b", mid = "gray", high = "#a70000", midpoint = 0, oob = scales::squish)
    scale_fill_custom = ggplot2::scale_fill_gradientn(limits = c(-6,6), breaks=seq(-6, 6, by = 3), colours = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100), oob = scales::squish)
    
    # scale_fill_custom = ggplot2::scale_fill_gradient2(limits = c(-6,6), breaks=seq(-6, 6, by = 2), low = "#00006b", mid = "white", high = "#a70000", midpoint = 0, oob = scales::squish)
    
  }
  
  
  
  new.seurat.object = CreateSeuratObject(counts = corrected_mix_mat_transposed, assay = "Spatial")
  
  new.seurat.object@images$image = new(
    
    Class = 'VisiumV1'
    
    ,assay = "spatial"
    
    ,key = "image_"
    
    ,coordinates = mapping_file
    
    ,image = visium@images[[image_id]]@image
    
    ,scale.factors = visium@images[[image_id]]@scale.factors
    
    ,spot.radius = visium@images[[image_id]]@spot.radius
    
  )
  
  # pt_size = round(max(max(mapping_file$row), max(mapping_file$col))/120, 2) ### does not work
  
  
  p = SpatialFeaturePlot(new.seurat.object
                         
                         , features = tc_to_plot
                         ,slot = "counts"
                         , stroke = NA
                         
                         , image.alpha = 0.2
                         
                         # , pt.size.factor = 1.8
                         
                         # , pt.size.factor = 1.4  ### for OC1 zoomed
                         
                         # , pt.size.factor = pt_size
                         
                         , alpha = 1
                         
                         ,interactive = FALSE)&scale_fill_custom&NoAxes()
  
  # ,interactive = FALSE)&NoAxes()
  
  p
  
}



############################################

### iterate through all

images_to_plot = c("BC1", "BC2", "BC3", "BC4", "BC5", "CC1", "CRC1", "CRC2", "CRC3", "GB1", "GB2", "IDC1", "LC1", "LC2", "M1", "OC1", "OC2", "OC3", "PC1", "PC2", "PC3")



plot_list_tcs = lapply(images_to_plot, function(curr_img){
  
  print(paste0("### start ", curr_img))
  
  plot_list_tcs_curr_tc = lapply(c("text", 'raw',tcs_to_plot), function(curr_tc){
    
    print(paste0("### ", curr_tc))
    
    p = NULL
    
    if(curr_tc =='text') {
      
      p = ggplot() + geom_text(aes(label = I(curr_img), y=0.5, x=0.2), hjust = 0) + xlim(0, 1) + ylim(0,1) + theme_void()
      
    } else {
      
      if(curr_tc=='raw'){
        
        p = plot_raw_image(curr_img, CNA_corrected_mix_mat = CNA_corrected_mix_mat, visium = visium)
        
        p = p + theme(plot.margin=margin(l=-0.5,unit="in"))
        
      } else {
        
        p = plot_tc_image(curr_img, tc_to_plot=curr_tc,CNA_corrected_mix_mat = CNA_corrected_mix_mat, visium = visium)
        
        p = p + theme(plot.margin=margin(l=-0.5,unit="in"))
        
      }
      
    }
    
    p
    
  })
  
  plot_list_tcs_curr_tc
  
})

plot_list_tcs = unlist(plot_list_tcs, recursive = F)



pdf(paste0('TC_TCGA_100EV_25CNA_TCs.pdf'), width = 3*length(tcs_to_plot), height = 2*length(images_to_plot))

plot_grid(plotlist = plot_list_tcs, nrow=length(images_to_plot), rel_widths = c(0.6,rep(1,length(tcs_to_plot)+2)), scale = 1)

dev.off()




