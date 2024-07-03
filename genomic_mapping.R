################## Genomic mapping in CNA-TCs #######################
### genomic ordered plot of tcs

library(data.table)
library(ggplot2)
library(dplyr)
library(glue)


## set the working directory in R to the directory where the currently active script is located.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))   ### only in RStudio


WEIGHT_LIMIT <- 20 ### cutoff for plot
PLOT_RECT_ALPHA <- 0.07

OUTPUT_DIR = getwd()

plot_dir = file.path(OUTPUT_DIR, "genomic_plots")
dir.create(plot_dir)

##### TACNA_CNA_TCs in original flipped file: cna_ica_flipped_df #######

############################################
### read in

gene_anno = fread("/home/iro/CCLE/Entrezid_mapping_using_org_Hs_eg_db_17012024.txt")
gene_anno$mapped_entrez_v1 = as.character(gene_anno$mapped_entrez_v1)

columns_to_plot = fread('/home/iro/GSEA_results/cna_ica_flipped_df.tsv')

# select the TCs to plot 
columns_to_plot = columns_to_plot[,]


############################################
### preprocessing

# columns_to_plot = data.table(t(columns_to_plot), keep.rownames = T)
colnames(columns_to_plot)[1] = "V1"
columns_to_plot$V1 = as.character(columns_to_plot$V1)

columns_to_plot[,2:ncol(columns_to_plot)] = pmax(columns_to_plot[,2:ncol(columns_to_plot)], -(WEIGHT_LIMIT-0.5))
columns_to_plot[,2:ncol(columns_to_plot)] = pmin(columns_to_plot[,2:ncol(columns_to_plot)], (WEIGHT_LIMIT-0.5))


### edit to to gene_id|chr|base
mapping_df = gene_anno %>%
  filter(chromosome_no %in% c(1:22,'X','Y')) %>%
  dplyr::mutate("bp_middle" = (BP_startloc + (BP_endloc-BP_startloc)/2) ) %>%
  dplyr::select(mapped_entrez_v1, chromosome_no, bp_middle)

# table(mapping_df$chr)
colnames(mapping_df) = c("gene_id", "chr_mapping","bp_mapping")
mapping_df = data.table(mapping_df)

plot_file <- mapping_df[columns_to_plot, on="gene_id==V1"] 


############################################
### edit for plotting input

plot_file[chr_mapping=="", chr_mapping := NA]
plot_file[chr_mapping=="X", chr_mapping := "23"]
plot_file[chr_mapping=="Y", chr_mapping := "24"]
plot_file[, chr := as.numeric(chr_mapping)]
plot_file = plot_file %>% tidyr::drop_na(chr)

plot_file[order(chr,bp_mapping),plot_order := .I]

chromosome_annotation <- plot_file[,.(median(.SD[,plot_order]) %>% round(),min(plot_order),max(plot_order)), by=chr]
setnames(chromosome_annotation, c("chr", "middle","start","end"))
chromosome_annotation[,chr_text := as.character(chr)]
chromosome_annotation[chr==23,chr_text := "X"]
chromosome_annotation[chr==24,chr_text := "Y"]


############################################
### actual plotting


plot_ic_genome_ordered_weights <- function(component, plot_file, output_dir, main=""){
  
  print(component)
  comp_plot_df = plot_file[!is.na(chr)][order(plot_order)]
  
  p <- ggplot(comp_plot_df, aes_string("plot_order", "`{component}`"%>%glue)) + 
    annotate("rect", xmin = chromosome_annotation[chr==1,start], xmax = chromosome_annotation[chr==1,end], ymin = -WEIGHT_LIMIT, ymax = WEIGHT_LIMIT, alpha = PLOT_RECT_ALPHA) +
    annotate("rect", xmin = chromosome_annotation[chr==3,start], xmax = chromosome_annotation[chr==3,end], ymin = -WEIGHT_LIMIT, ymax = WEIGHT_LIMIT, alpha = PLOT_RECT_ALPHA) +
    annotate("rect", xmin = chromosome_annotation[chr==5,start], xmax = chromosome_annotation[chr==5,end], ymin = -WEIGHT_LIMIT, ymax = WEIGHT_LIMIT, alpha = PLOT_RECT_ALPHA) +
    annotate("rect", xmin = chromosome_annotation[chr==7,start], xmax = chromosome_annotation[chr==7,end], ymin = -WEIGHT_LIMIT, ymax = WEIGHT_LIMIT, alpha = PLOT_RECT_ALPHA) +
    annotate("rect", xmin = chromosome_annotation[chr==9,start], xmax = chromosome_annotation[chr==9,end], ymin = -WEIGHT_LIMIT, ymax = WEIGHT_LIMIT, alpha = PLOT_RECT_ALPHA) +
    annotate("rect", xmin = chromosome_annotation[chr==11,start], xmax = chromosome_annotation[chr==11,end], ymin = -WEIGHT_LIMIT, ymax = WEIGHT_LIMIT, alpha = PLOT_RECT_ALPHA) +
    annotate("rect", xmin = chromosome_annotation[chr==13,start], xmax = chromosome_annotation[chr==13,end], ymin = -WEIGHT_LIMIT, ymax = WEIGHT_LIMIT, alpha = PLOT_RECT_ALPHA) +
    annotate("rect", xmin = chromosome_annotation[chr==15,start], xmax = chromosome_annotation[chr==15,end], ymin = -WEIGHT_LIMIT, ymax = WEIGHT_LIMIT, alpha = PLOT_RECT_ALPHA) +
    annotate("rect", xmin = chromosome_annotation[chr==17,start], xmax = chromosome_annotation[chr==17,end], ymin = -WEIGHT_LIMIT, ymax = WEIGHT_LIMIT, alpha = PLOT_RECT_ALPHA) +
    annotate("rect", xmin = chromosome_annotation[chr==19,start], xmax = chromosome_annotation[chr==19,end], ymin = -WEIGHT_LIMIT, ymax = WEIGHT_LIMIT, alpha = PLOT_RECT_ALPHA) +
    annotate("rect", xmin = chromosome_annotation[chr==21,start], xmax = chromosome_annotation[chr==21,end], ymin = -WEIGHT_LIMIT, ymax = WEIGHT_LIMIT, alpha = PLOT_RECT_ALPHA) +
    annotate("rect", xmin = chromosome_annotation[chr==23,start], xmax = chromosome_annotation[chr==23,end], ymin = -WEIGHT_LIMIT, ymax = WEIGHT_LIMIT, alpha = PLOT_RECT_ALPHA) +
    geom_point(aes(alpha = 0.2), pch='.', size=2) +
    scale_color_manual(values = c("black")) +
    # geom_point(aes(col=chr%%2==0, alpha = 0.1), size=0.00001) +
    # scale_color_manual(values = c("#307bc5", "#cfa700")) +
    ylab(" IC gene weights") +
    xlab("genomic mapping") +
    ggtitle(paste0(main, component)) +
    scale_x_continuous(expand = c(0, 0), breaks = chromosome_annotation[,middle], labels = chromosome_annotation[,chr_text]) + 
    scale_y_continuous(expand = c(0, 0), limits = c(-WEIGHT_LIMIT, +WEIGHT_LIMIT)) +
    theme_bw() +
    theme(axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size = 10),
          plot.title = element_text(size = 10),
          panel.border = element_rect(colour = alpha('black', 0.2), fill = NA),
    ) + guides(color=FALSE, alpha=FALSE)
  
  underscore_name <- gsub(" ", "_", component)
  ggsave(file.path(output_dir, "{underscore_name}.png") %>% glue, p, dpi = 300, height = 7, width = 30, units = "cm" )
}


### plot all components
for (component in colnames(columns_to_plot[,-1])){
  plot_ic_genome_ordered_weights(component, plot_file=plot_file, output_dir=plot_dir)
}
