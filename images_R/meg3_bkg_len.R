
library(ggplot2)
library(cowplot)  # To have ggplots side-by-side: plot_grid()

###

DATA_DIR <- "../data/"
OUT_DIR <- "../images/"

theme_set(theme_bw(base_size = 19))  # increase the font size: https://stackoverflow.com/a/11955412/310453

###

get_ggplot_from_fn <- function(bed_fn, name)
{
  # bed_fn <- 'MEG3_peaks_hg38.bed'
  # name <- 'MEG3 peaks'
  bed_df <- read.delim(paste0(DATA_DIR, bed_fn), as.is=TRUE, header=FALSE)
  colnames(bed_df)[1:3] <- c('chr','start','end')
  bed_df$len <- bed_df$end - bed_df$start
  nrow(bed_df)
  
  title <- sprintf("%s (n=%d)", name, nrow(bed_df))
  subT <- sprintf("median length = %.1f bp, longest peak = %d bp",
                  median(bed_df$len), max(bed_df$len))
  ggplot(bed_df) +
    aes(x = len) +
    geom_histogram(binwidth = 30, col = "black", fill = "white") +
    xlim(0, 2000) +
    labs(x = "Peak length (bp)", y = "Number of peaks") +
    ggtitle(title, subtitle = subT)
}

meg3_gg <- get_ggplot_from_fn('MEG3_peaks_hg38.bed', 'MEG3 peaks')
bkg_gg <- get_ggplot_from_fn('Background_hg38.bed', 'Control DNA regions')

plot_grid(meg3_gg, bkg_gg, ncol = 1, labels = 'auto')
ggsave('meg3_bkg_len.pdf', path=OUT_DIR, width = 8, height = 8)
