
############################################
### create colocalization corr between TCs

# setwd(dirname(rstudioapi::getSourceEditorContext()$path))   ### only in RStudio

source("https://raw.githubusercontent.com/loipf/coding_snippets/master/R/small_functions.R")
source("https://raw.githubusercontent.com/loipf/coding_snippets/master/R/save_load_dataframe.R")

options(bitmapType='cairo')

#library(Seurat)
library(BiocParallel)
library(ks)
library(dplyr)
library(SeuratObject)
library(parallel)
library(parallelly)
no_cores <- availableCores()
# ### adapted from:
# https://github.com/sameelab/STANN/blob/main/downstream_analysis/run_kde_estimates.R
# https://github.com/sameelab/STANN/blob/main/downstream_analysis/run_pearson_correlation.py



############################################
### read in 
setwd("/scratch/hb-bioinfo-fehrmann/projects_loipf/coloc_iro/")
corrected_mix_mat = data.frame(fread('corrected_mix_mat.tsv', stringsAsFactors = F), row.names = 1, check.names = F)

visium = readRDS("ALL_VISIUM_FILES_MERGED.rds")

# BBPARAM_multicore = MulticoreParam(10)  ### local
BBPARAM_multicore = BiocParallel::SnowParam(24, type = "FORK") ### cluster

# OUTPUT_DIR = file.path("/scratch/hb-bioinfo-fehrmann/projects_loipf/coloc_iro/")
OUTPUT_DIR = file.path("/scratch/hb-bioinfo-fehrmann/Projects_from_Arko/Example_directory/")


############################################
### functions

create_mapping_mm = function(image_id, corrected_mix_mat = corrected_mix_mat_01, visium = visium){
  mapping_file <- visium@images[[image_id]]@coordinates
  
  common_samples = intersect(rownames(mapping_file),rownames(corrected_mix_mat) )
  mapping_file = mapping_file[common_samples,]
  corrected_mix_mat_current = corrected_mix_mat[common_samples,]
  mapping_mm_df = cbind(mapping_file, corrected_mix_mat_current)
  mapping_mm_df
}

matrix_quantile_cutoff = function(x, q=0.5) {
  quantile_cutoff = quantile(x, q)
  x>quantile_cutoff
}



############################################
### calculate kde estimate

images_to_plot = grep("EC", names(visium@images), value = T, invert = T)  ### EC is inhouse
#images_to_plot = c("BC1", "BC2", "BC3", "BC4", "BC5", "CC1", "CRC1", "CRC2", "CRC3", "GB1", "GB2", "IDC1", "LC1", "LC2", "M1", "OC1", "OC2", "OC3", "PC1", "PC2", "PC3")

tcs_to_plot = colnames(corrected_mix_mat)

### which side
curr_side = 'neg2'
# curr_side = 'pos2'


# job::job({
  dir.create(file.path(OUTPUT_DIR, paste0("ks_plots",curr_side)), showWarnings = F)
  img_tc_kernels_positive = sapply(images_to_plot, function(curr_img) {
    print(curr_img)
    curr_mapping_df_mm = create_mapping_mm(curr_img, corrected_mix_mat = corrected_mix_mat, visium = visium)
    # curr_mapping_df_mm$cnaburden_per_image = range01(curr_mapping_df_mm$cnaburden_01)

    OUTPUT_DIR_IMG = file.path(OUTPUT_DIR, paste0("ks_plots",curr_side), curr_img)
    dir.create(OUTPUT_DIR_IMG, showWarnings = F)

    ### set origin to 0
    curr_mapping_df_mm$row = curr_mapping_df_mm$row - min(curr_mapping_df_mm$row)
    curr_mapping_df_mm$col = curr_mapping_df_mm$col - min(curr_mapping_df_mm$col)

    curr_mapping_limits_row = c(min(curr_mapping_df_mm$row), max(curr_mapping_df_mm$row))
    curr_mapping_limits_col = c(min(curr_mapping_df_mm$col), max(curr_mapping_df_mm$col))
    x_kde_grid = matrix(NA, nrow = curr_mapping_limits_row[2]-curr_mapping_limits_row[1], ncol = curr_mapping_limits_col[2]-curr_mapping_limits_col[1])

    # tc_kde_list = list()
    # for(i in 1:length(tcs_to_plot))
    # {
    #   curr_tc = tcs_to_plot[i]
    #   if(curr_tc =="cnaburden_01") {
    #     x = curr_mapping_df_mm[curr_mapping_df_mm[['cnaburden_01']] > CNABURDEN_CUTOFF ,]
    #     # x = curr_mapping_df_mm[curr_mapping_df_mm[['cnaburden_01']] < CNABURDEN_CUTOFF ,]
    #   } else {
    #     if(curr_side =='pos2'){
    #       mm_cutoff = 2
    #       x = curr_mapping_df_mm[curr_mapping_df_mm[[curr_tc]] > mm_cutoff ,]  ### pos only
    #     } else {
    #       if (curr_side =='neg2'){
    #         mm_cutoff = -2
    #         x = curr_mapping_df_mm[curr_mapping_df_mm[[curr_tc]] < mm_cutoff , ]  ### neg only
    #       }
    #     }
    #   }
    #   
    #   ### at least n spots must be there to identify pattern
    #   if(nrow(x)<50){
    #     # if(nrow(x)<30){
    #     x_kde = x_kde_grid
    #     tc_kde_list[[i]] = x_kde
    #     print(i)
    #   } else {
    #     
    #     ### bandwith optimisation start - avoids bandwith close to 0
    #     h_start <- matrix(0, nrow = 2, ncol = 2)
    #     diag(h_start) <- 9
    #     
    #     ### Hlscv.diag uses all threads - no idea how to limit
    #     x_hpi = ks::Hlscv.diag(x=x[,c("row","col")], Hstart=h_start, bgridsize = dim(x_kde_grid))
    #     diag(x_hpi) = pmax(diag(x_hpi),3)  ### rare edge cases: sometimes bandwith optimisation is close to 0
    #     
    #     ### kernel fit
    #     x_kde = ks::kde(x=x[,c("row","col")], w=abs(x[,curr_tc]), H=x_hpi, binned=T, bgridsize = dim(x_kde_grid), density = TRUE, xmin = c(0,0), xmax = dim(x_kde_grid))
    #     
    #     ### plot output
    #     png(file.path(OUTPUT_DIR_IMG, paste0(curr_tc,".png")), width=1050, height = 400)
    #     par(mfrow=c(1,3))
    #     plot(x$row, x$col, ann=FALSE,axes=F,frame.plot=T, ylim=curr_mapping_limits_col, xlim=curr_mapping_limits_row)
    #     image(x_kde$estimate, useRaster = T, axes = FALSE)
    #     
    #     x_kde_mask = matrix_quantile_cutoff(x_kde$estimate, 0.75)
    #     x_kde_plot = x_kde$estimate
    #     x_kde_plot[!x_kde_mask] = NA
    #     image(x_kde_plot)
    #     dev.off()
    #     
    #   tc_kde_list[[i]] = x_kde$estimate
    #   print(i)
    #   }
    # }
    
    
    cl = makeCluster(5, type = "FORK")
    
    tc_kde_list = parLapply(cl,tcs_to_plot, function(curr_tc) {
      if(curr_tc =="cnaburden_01") {
        x = curr_mapping_df_mm[curr_mapping_df_mm[['cnaburden_01']] > CNABURDEN_CUTOFF ,]
        # x = curr_mapping_df_mm[curr_mapping_df_mm[['cnaburden_01']] < CNABURDEN_CUTOFF ,]
      } else {
        if(curr_side =='pos2'){
          mm_cutoff = 2
          x = curr_mapping_df_mm[curr_mapping_df_mm[[curr_tc]] > mm_cutoff ,]  ### pos only
        } else {
          if (curr_side =='neg2'){
            mm_cutoff = -2
            x = curr_mapping_df_mm[curr_mapping_df_mm[[curr_tc]] < mm_cutoff , ]  ### neg only
          }
        }
      }
      
      ### at least n spots must be there to identify pattern
      if(nrow(x)<50){
        # if(nrow(x)<30){
        x_kde = x_kde_grid
        return(x_kde)
      } else {
        
        ### bandwith optimisation start - avoids bandwith close to 0
        h_start <- matrix(0, nrow = 2, ncol = 2)
        diag(h_start) <- 9
        
        ### Hlscv.diag uses all threads - no idea how to limit
        x_hpi = ks::Hlscv.diag(x=x[,c("row","col")], Hstart=h_start, bgridsize = dim(x_kde_grid))
        diag(x_hpi) = pmax(diag(x_hpi),3)  ### rare edge cases: sometimes bandwith optimisation is close to 0
        
        ### kernel fit
        x_kde = ks::kde(x=x[,c("row","col")], w=abs(x[,curr_tc]), H=x_hpi, binned=T, bgridsize = dim(x_kde_grid), density = TRUE, xmin = c(0,0), xmax = dim(x_kde_grid))
        
        ### plot output
        png(file.path(OUTPUT_DIR_IMG, paste0(curr_tc,".png")), width=1050, height = 400)
        par(mfrow=c(1,3))
        plot(x$row, x$col, ann=FALSE,axes=F,frame.plot=T, ylim=curr_mapping_limits_col, xlim=curr_mapping_limits_row)
        image(x_kde$estimate, useRaster = T, axes = FALSE)
        
        x_kde_mask = matrix_quantile_cutoff(x_kde$estimate, 0.75)
        x_kde_plot = x_kde$estimate
        x_kde_plot[!x_kde_mask] = NA
        image(x_kde_plot)
        dev.off()
        
        return(x_kde$estimate)
      }
      
    })
    
    stopCluster(cl)
    
    
    
    
    # tc_kde_list = lapply(tcs_to_plot, function(curr_tc) {
    # tc_kde_list = bplapply(tcs_to_plot, function(curr_tc) {
    #   if(curr_tc =="cnaburden_01") {
    #     x = curr_mapping_df_mm[curr_mapping_df_mm[['cnaburden_01']] > CNABURDEN_CUTOFF ,]
    #     # x = curr_mapping_df_mm[curr_mapping_df_mm[['cnaburden_01']] < CNABURDEN_CUTOFF ,]
    #   } else {
    #     if(curr_side =='pos2'){
    #       mm_cutoff = 2
    #       x = curr_mapping_df_mm[curr_mapping_df_mm[[curr_tc]] > mm_cutoff ,]  ### pos only
    #     } else {
    #       if (curr_side =='neg2'){
    #         mm_cutoff = -2
    #         x = curr_mapping_df_mm[curr_mapping_df_mm[[curr_tc]] < mm_cutoff , ]  ### neg only
    #       }
    #     }
    #   }
    # 
    #   ### at least n spots must be there to identify pattern
    #   if(nrow(x)<50){
    #   # if(nrow(x)<30){
    #     x_kde = x_kde_grid
    #     return(x_kde)
    #   } else {
    # 
    #     ### bandwith optimisation start - avoids bandwith close to 0
    #     h_start <- matrix(0, nrow = 2, ncol = 2)
    #     diag(h_start) <- 9
    # 
    #     ### Hlscv.diag uses all threads - no idea how to limit
    #     x_hpi = ks::Hlscv.diag(x=x[,c("row","col")], Hstart=h_start, bgridsize = dim(x_kde_grid))
    #     diag(x_hpi) = pmax(diag(x_hpi),3)  ### rare edge cases: sometimes bandwith optimisation is close to 0
    # 
    #     ### kernel fit
    #     x_kde = ks::kde(x=x[,c("row","col")], w=abs(x[,curr_tc]), H=x_hpi, binned=T, bgridsize = dim(x_kde_grid), density = TRUE, xmin = c(0,0), xmax = dim(x_kde_grid))
    # 
    #     ### plot output
    #     png(file.path(OUTPUT_DIR_IMG, paste0(curr_tc,".png")), width=1050, height = 400)
    #     par(mfrow=c(1,3))
    #     plot(x$row, x$col, ann=FALSE,axes=F,frame.plot=T, ylim=curr_mapping_limits_col, xlim=curr_mapping_limits_row)
    #     image(x_kde$estimate, useRaster = T, axes = FALSE)
    # 
    #     x_kde_mask = matrix_quantile_cutoff(x_kde$estimate, 0.75)
    #     x_kde_plot = x_kde$estimate
    #     x_kde_plot[!x_kde_mask] = NA
    #     image(x_kde_plot)
    #     dev.off()
    # 
    #     return(x_kde$estimate)
    #   }
    # 
    # },BPPARAM = BBPARAM_multicore)
    # # })

    names(tc_kde_list) = tcs_to_plot
    tc_kde_list

  }, USE.NAMES = T, simplify=F)

  saveRDS(img_tc_kernels_positive, file.path(OUTPUT_DIR, paste0("spatial_kernels_tc_",curr_side,"_obj.rds")))

# })


  
  
  
# ############################################
# ### calculate colocalization correlation
# 
# print_time('read in spatial_kernal obj rds')
# 
# ### calc corr all pos|neg possibilties ~4 hours
# # job::job({
#   img_tc_kernels = list('pos'=readRDS(file.path(OUTPUT_DIR, "spatial_kernels_tc_pos2_obj.rds")),
#                         'neg'=readRDS(file.path(OUTPUT_DIR, "spatial_kernels_tc_neg2_obj.rds")) )
#   curr_img=names(img_tc_kernels$pos)[1] ### just need one example
# 
#   par_grid = expand.grid('image_id'=names(img_tc_kernels$pos),
#                          'tc_1' = names(img_tc_kernels$pos[[curr_img]]),'tc_2' = names(img_tc_kernels$pos[[curr_img]]),
#                          'tc1_direction' = c('pos','neg'), 'tc2_direction' = c('pos','neg'), stringsAsFactors=F)
#   # par_grid = subset(par_grid, tc_1 =='cnaburden_01')
# 
# print_time('start corr calc')
# 
#   corr_list = lapply(1:nrow(par_grid), function(curr_row) {
#   # corr_list = bplapply(1:nrow(par_grid), function(curr_row) {
#     # print(curr_row)
#     curr_par_grid = par_grid[curr_row,]
# 
#     if(curr_row %%10000 == 0) {
#          print_time(paste0(curr_row,' / ',nrow(par_grid)) )
#     }
# 
#     tc_1 = img_tc_kernels[[curr_par_grid$tc1_direction]][[curr_par_grid$image_id]][[curr_par_grid$tc_1]]
#     tc_2 = img_tc_kernels[[curr_par_grid$tc2_direction]][[curr_par_grid$image_id]][[curr_par_grid$tc_2]]
# 
#     if(all(is.na(tc_1)) | all(is.na(tc_2))) {
#       curr_par_grid$pearson_corr = NA
#       curr_par_grid$pearson_corr_pvalue = NA
#       return(curr_par_grid)
#     }
# 
#     tc_1_mask = matrix_quantile_cutoff(tc_1, 0.75)  ### 0.75 not like 0.5 like in original
#     tc_2_mask = matrix_quantile_cutoff(tc_2, 0.75)
#     kde_mask_common = tc_1_mask | tc_2_mask  ### union not intersect like in original
#     corr_obj = cor.test(tc_1[kde_mask_common], tc_2[kde_mask_common], method="pearson")
# 
#     curr_par_grid$pearson_corr = corr_obj$estimate
#     curr_par_grid$pearson_corr_pvalue = corr_obj$p.value
#     return(curr_par_grid)
#   # }, BPPARAM = BBPARAM_multicore)
#   })
#   corr_df = do.call('rbind',corr_list)
#   save_dataframe(corr_df, file.path(OUTPUT_DIR,paste0('spatial_kernels_correlation.tsv')))
# # })
# 
# 
# print_time('end corr calc')
# 
# 
# ############################################
# ### calc combi coloc-exclusive
# 
# spatial_kernel_corr_df = load_dataframe(file.path(OUTPUT_DIR,paste0('spatial_kernels_correlation.tsv')))
#     
# spatial_kernel_corr_df_sig = spatial_kernel_corr_df
# spatial_kernel_corr_df_sig$pearson_corr[spatial_kernel_corr_df_sig$pearson_corr_pvalue>0.05] = NA
# spatial_kernel_corr_df_sig$pearson_corr[is.na(spatial_kernel_corr_df_sig$pearson_corr)] = 0
# # spatial_kernel_corr_df_sig = subset(spatial_kernel_corr_df_sig, tc_1 %in% cna_ic_list & tc_2 %in% immune_ic_list)
# spatial_kernel_corr_df_sig$unique_tc_id = paste0(spatial_kernel_corr_df_sig$image_id,'_',spatial_kernel_corr_df_sig$tc_1,'_',spatial_kernel_corr_df_sig$tc_2)
# 
# print_time('start calc coloc score4')
# 
# ## coloc_score = (pos_pos - pos_neg) - (neg_pos - neg_neg)
# for(curr_img in unique(spatial_kernel_corr_df_sig$image_id)) {
#   
#   print_time(curr_img)
#   curr_spatial_kernel_corr_df_sig = subset(spatial_kernel_corr_df_sig, image_id == curr_img)
#   
#   coloc_score_list = bplapply(unique(curr_spatial_kernel_corr_df_sig$unique_tc_id), function(curr_id) {
#     curr_df = subset(curr_spatial_kernel_corr_df_sig, unique_tc_id==curr_id)
#     coloc_pos_pos = subset(curr_df, tc1_direction=="pos" & tc2_direction=="pos")$pearson_corr
#     coloc_pos_neg = subset(curr_df, tc1_direction=="pos" & tc2_direction=="neg")$pearson_corr
#     coloc_neg_pos = subset(curr_df, tc1_direction=="neg" & tc2_direction=="pos")$pearson_corr
#     coloc_neg_neg = subset(curr_df, tc1_direction=="neg" & tc2_direction=="neg")$pearson_corr
#     coloc_score = (coloc_pos_pos - coloc_pos_neg) - (coloc_neg_pos - coloc_neg_neg)
#     
#     ### division by number of available non-NA values
#     number_non_NAs = sum(c(coloc_pos_pos,coloc_pos_neg, coloc_neg_pos,coloc_neg_neg)!=0)
#     # coloc_score = coloc_score/number_non_NAs
#     
#     data.frame("tc_id"=curr_id, "coloc_score" = coloc_score,'number_non_NAs'=number_non_NAs)
#     # })
#   }, BPPARAM = BBPARAM_multicore)
#   # })
#   coloc_score_df = do.call("rbind",coloc_score_list)
#   coloc_score_df = coloc_score_df %>% tidyr::separate(tc_id, into = c("image_id","tc_1","tc_2"), sep = '_')
#   coloc_score_df$tc_id = gsub("consensus independent component ","TC",paste0(coloc_score_df$tc_1,'_',coloc_score_df$tc_2))
#   save_dataframe(coloc_score_df, file.path(OUTPUT_DIR, paste0('spatial_kernels_coloc_score4_cutoff2_',curr_img,'.tsv')))
# }
# 
# 
# combi_coloc_score_df_list = lapply(unique(spatial_kernel_corr_df_sig$image_id), function(curr_img){
#   curr_df = load_dataframe(file.path(OUTPUT_DIR,paste0('spatial_kernels_coloc_score4_cutoff2_',curr_img,'.tsv')))
#   curr_df
# })
# combi_coloc_score_df = do.call('rbind', combi_coloc_score_df_list)
# save_dataframe(combi_coloc_score_df, file.path(OUTPUT_DIR, paste0('spatial_kernels_coloc_score4_cutoff2_all.tsv')))





# ### run coloc4score calculation all at once - memory issues 
# # coloc_score_list = lapply(unique(spatial_kernel_corr_df_sig$unique_tc_id), function(curr_id) {
# # job::job({
# 
#   coloc_score_list = bplapply(unique(spatial_kernel_corr_df_sig$unique_tc_id), function(curr_id) {
#   # coloc_score_list = lapply(unique(spatial_kernel_corr_df_sig$unique_tc_id), function(curr_id) {
#     curr_df = subset(spatial_kernel_corr_df_sig, unique_tc_id==curr_id)
#     coloc_pos_pos = subset(curr_df, tc1_direction=="pos" & tc2_direction=="pos")$pearson_corr
#     coloc_pos_neg = subset(curr_df, tc1_direction=="pos" & tc2_direction=="neg")$pearson_corr
#     coloc_neg_pos = subset(curr_df, tc1_direction=="neg" & tc2_direction=="pos")$pearson_corr
#     coloc_neg_neg = subset(curr_df, tc1_direction=="neg" & tc2_direction=="neg")$pearson_corr
#     coloc_score = (coloc_pos_pos - coloc_pos_neg) - (coloc_neg_pos - coloc_neg_neg)
# 
#     ### division by number of available non-NA values
#     number_non_NAs = sum(c(coloc_pos_pos,coloc_pos_neg, coloc_neg_pos,coloc_neg_neg)!=0)
#     # coloc_score = coloc_score/number_non_NAs
# 
#     data.frame("tc_id"=curr_id, "coloc_score" = coloc_score,'number_non_NAs'=number_non_NAs)
#     # })
#   }, BPPARAM = BBPARAM_multicore)
# # })
# coloc_score_df = do.call("rbind",coloc_score_list)
# coloc_score_df = coloc_score_df %>% tidyr::separate(tc_id, into = c("image_id","tc_1","tc_2"), sep = '_')
# coloc_score_df$tc_id = gsub("consensus independent component ","TC",paste0(coloc_score_df$tc_1,'_',coloc_score_df$tc_2))
# save_dataframe(coloc_score_df, 'spatial_kernels_coloc_score4_cutoff2_parallel32.tsv')
# save_dataframe(coloc_score_df, paste0('spatial_kernels_coloc_score4_cutoff2_',curr_img,'.tsv'))
# 
# 




print_time('finished')


