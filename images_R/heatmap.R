
library(ggplot2)
library(plyr)   # revalue()
library(dplyr)

library(ComplexHeatmap)
library(circlize)       # colorRamp2()

source("lib.R")   # get_log_pvalue_mtx(), find_uni_tts()


###

DATA_DIR <- "../data/"
OUT_DIR <- "../images/"

###

peak_init_df <- read.table(paste0(DATA_DIR, "pvalue_chop_seq.txt.gz"), header=TRUE, row.names=1)
bkg_init_df <- read.table(paste0(DATA_DIR, "pvalue_background.txt.gz"), header=TRUE, row.names=1)
capture_init_df <- read.table(paste0(DATA_DIR, "pvalue_capture_seq.txt.gz"), header=TRUE, row.names=1)

# Are all p-value valid? (must be TRUE)
all(0 <= peak_init_df & peak_init_df <= 1)
all(0 <= bkg_init_df & bkg_init_df <= 1)
all(0 <= capture_init_df & capture_init_df <= 1)

peak_log_p <- get_log_pvalue_mtx(peak_init_df)
bkg_log_p <- get_log_pvalue_mtx(bkg_init_df)
capture_log_p <- get_log_pvalue_mtx(capture_init_df)

peak_gg_df <- find_uni_tts(peak_log_p)
bkg_gg_df <- find_uni_tts(bkg_log_p)
capture_gg_df <- find_uni_tts(capture_log_p)


get_heatmap <- function(mtx, ...)
{
  return(
    Heatmap(mtx,
            row_title_gp = gpar(fontsize = 10),
            column_title_gp = gpar(fontsize = 10),
            show_column_names = FALSE,
            show_row_names = FALSE,
            show_column_dend = FALSE,
            show_row_dend = FALSE,
            col = colorRamp2(c(0, 3, 6, 10, 100), c("blue", "white", "yellow", "orange", "red")),
            heatmap_legend_param = list(legend_direction = "horizontal", title_position = "lefttop", legend_width = unit(5, "cm")),
            ...)
  )
}

get_three_heatmaps <- function(all.data, clust_v)
{
  # The cluster names are the number of the DNA regions in each cluster
  # http://www.cookbook-r.com/Manipulating_data/Renaming_levels_of_a_factor/
  clust_f <- revalue(factor(clust_v), c(
    '1_uni_tts'=as.character(sum(clust_v == '1_uni_tts')),
    '2_other'=as.character(sum(clust_v == '2_other'))))
  
  # Make the the main heatmap: trxs_ht
  trxs_mtx <- as.matrix(all.data[,grep('^ENST', colnames(all.data))])
  trxs_ht <- get_heatmap(trxs_mtx,
                         name = "-log10(adj_pvalue)",
                         column_title = sprintf("%d expressed RNAs", ncol(trxs_mtx)),
                         split = clust_f,
                         # combined_name_fun = NULL,   # Do not show cluster names
                         gap = unit(2, "mm"))
  
  # rand2_ht
  rand2_mtx <- as.matrix(all.data[,grep('^rand_', colnames(all.data))])
  rand2_ht <- get_heatmap(rand2_mtx,
                          column_title = sprintf("%d shuffled MEG3", ncol(rand2_mtx)),
                          show_heatmap_legend = FALSE)
  
  # meg3_ht
  meg3_mtx <- as.matrix(all.data[,grep('^MEG3$', colnames(all.data))])
  meg3_ht <- get_heatmap(meg3_mtx,
                         column_title = "MEG3",
                         show_heatmap_legend = FALSE,
                         width = unit(5, "mm"))
  
  return(trxs_ht + meg3_ht + rand2_ht)
}

save_three_heatmaps <- function(data, clust_v, title, png_fn)
{
  png(png_fn, width = 1000, height = 1000, units = "px", res = 200)
  draw(get_three_heatmaps(data, clust_v),
       heatmap_legend_side = "bottom",
       gap = unit(2, "mm"),
       row_title = paste(nrow(data), title),
       column_title = "Query RNAs")
  dev.off()
}

save_three_heatmaps(peak_log_p, peak_gg_df$type,
                    "ChOP-seq peaks", paste0(OUT_DIR, 'heatmap_chop.png'))
save_three_heatmaps(bkg_log_p, bkg_gg_df$type,
                    "Control DNA regions", paste0(OUT_DIR, 'heatmap_bkg.png'))
save_three_heatmaps(capture_log_p, capture_gg_df$type,
                    "Shared Capture-seq peaks", paste0(OUT_DIR, 'heatmap_capture.png'))


###
# Testing on a subset of DNA regions
#quartz()

# data <- peak_log_p[sample(1:nrow(peak_log_p), 600),]
# clust_v <- find_uni_tts(data)$type
# title <- "ChOP-seq peaks"
# draw(get_three_heatmaps(data, clust_v),
#      heatmap_legend_side = "bottom",
#      gap = unit(2, "mm"),
#      row_title = paste(nrow(data), title),
#      column_title = "Query RNAs")
# 
# 
