
library(dplyr)
library(tidyr)  # gather()
library(ggplot2)
library(cowplot)   # To have ggplots side-by-side: plot_grid()

###

DATA_DIR <- "../data/"
OUT_DIR <- "../images/"
PVALUE_THR <- 0.01
FRACTION_THR <- 0.9

theme_set(theme_bw(base_size = 19))  # increase the font size: https://stackoverflow.com/a/11955412/310453

get_log_pvalue_mtx <- function(df)
{
  #df <- df[sample(1:nrow(df), 600),]
  # https://stackoverflow.com/a/32934725/310453
  pvalue_mtx <- df %>%
    as.matrix %>%
    as.vector %>% 
#    p.adjust(method='fdr') %>% 
    p.adjust(method="bonferroni") %>% 
    matrix(ncol=ncol(df))
  colnames(pvalue_mtx) <- colnames(df)
  rownames(pvalue_mtx) <- rownames(df)
  pvalue_mtx[pvalue_mtx < 1e-100] = 1e-100
  -log10(pvalue_mtx)
}

find_uni_tts <- function(mtx)
{
  trxs_mtx <- mtx[,grep('^ENST', colnames(mtx))]
  value_thr <- -log10(PVALUE_THR)
  vals <- apply(trxs_mtx, 1, function(x) sum(x > value_thr)/length(x))
  gg_df <- data.frame(dna = rownames(mtx),
                      q_fraction = vals,
                      type = ifelse(vals > FRACTION_THR, '1_uni_tts', '2_other'))
  return(gg_df)
}
