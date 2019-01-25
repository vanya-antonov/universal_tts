
library(dplyr)
library(tibble)    # rownames_to_column()
library(plyr)      # revalue() renames factor levels
library(ggplot2)
library(cowplot)  # To have ggplots side-by-side: plot_grid()

###

DATA_DIR <- "../data/"
OUT_DIR <- "../images/"

theme_set(theme_bw(base_size = 19))  # increase the font size: https://stackoverflow.com/a/11955412/310453

###

peak_df <- read.table(paste0(DATA_DIR, "pvalue_chop_seq.txt.gz"), header=TRUE, row.names=1)
bkg_df <- read.table(paste0(DATA_DIR, "pvalue_background.txt.gz"), header=TRUE, row.names=1)

adjust_pvalue_df <- function(df)
{
  # https://stackoverflow.com/a/32934725/310453
  pvalue_mtx <- as.matrix(df) %>%
    as.vector %>% 
#    p.adjust(method='fdr') %>% 
    p.adjust(method="bonferroni") %>% 
    matrix(ncol=ncol(df))
  colnames(pvalue_mtx) <- colnames(df)
  rownames(pvalue_mtx) <- rownames(df)
  as.data.frame(pvalue_mtx)
}

get_gg_df_from_pvalue_df <- function(df)
{
  adjust_pvalue_df(df) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column('RNA') %>%
    # https://stackoverflow.com/a/47591044/310453
    mutate(n_signif = rowSums(select(., -RNA) < 0.01)) %>%
    select(RNA, n_signif)
}

get_ggplot_from_gg_df <- function(gg_df, title)
{
  rna_types_v = c('1_meg3' = 'MEG3', '2_rna' = '153 expressed RNAs', '3_rand' = '153 shuffled MEG3')
  gg_df <- gg_df %>%
    mutate(RNA_type = case_when(RNA == 'MEG3' ~ '1_meg3',
                                startsWith(RNA, 'ENST') ~ '2_rna',
                                startsWith(RNA, 'rand') ~ '3_rand',
                                TRUE ~ '')) %>%
    mutate(RNA_type = factor(RNA_type)) %>%
    mutate(RNA_type = revalue(RNA_type, rna_types_v))
  
  subT = paste0('MEG3 interactions with adj p-value < 0.01 = ',
                filter(gg_df, RNA=='MEG3')$n_signif)
  ggplot(gg_df) +
    aes(x = RNA_type, y = n_signif) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    ylim(0, nrow(peak_df)) +
    xlab('Query RNA(s)') +
    ylab('Number of DNA regions\nwith adj p-value < 0.01') +
    ggtitle(title, subtitle = subT)
}

peak_gg_df <- get_gg_df_from_pvalue_df(peak_df)
bkg_gg_df <- get_gg_df_from_pvalue_df(bkg_df)

get_ggplot_from_gg_df(peak_gg_df, paste(nrow(peak_df), 'MEG3 ChOP-seq peaks'))
ggsave('num_signif_chop.pdf', path=OUT_DIR, width = 8, height = 5)

get_ggplot_from_gg_df(bkg_gg_df, paste(nrow(bkg_df), 'Control DNA regions'))
ggsave('num_signif_bkg.pdf', path=OUT_DIR, width = 8, height = 5)


# Statistics ----
total_peaks <- nrow(peak_df)
  
meg3_vs_peaks_n <- filter(peak_gg_df, RNA=='MEG3')$n_signif
meg3_vs_bkg_n <- filter(bkg_gg_df, RNA=='MEG3')$n_signif

# 3825 (56.3%)
sprintf('%d (%.1f%%)', meg3_vs_peaks_n, 100*meg3_vs_peaks_n/total_peaks)

# 617 (9.1%)
sprintf('%d (%.1f%%)', meg3_vs_bkg_n, 100*meg3_vs_bkg_n/total_peaks)

###

trxs_vs_peaks_gg_df <- filter(peak_gg_df, startsWith(RNA, 'ENST'))
trxs_vs_bkg_gg_df   <- filter(bkg_gg_df,  startsWith(RNA, 'ENST'))

# 65 (42%)
trxs_vs_peaks_n50 <- filter(trxs_vs_peaks_gg_df, n_signif > total_peaks/2) %>% nrow()
total_trxs <- nrow(trxs_vs_peaks_gg_df)
sprintf('%d (%.0f%%)', trxs_vs_peaks_n50, 100*trxs_vs_peaks_n50/total_trxs)

# 2890 (43%)
trxs_vs_peaks_median <- as.integer(median(trxs_vs_peaks_gg_df$n_signif))
sprintf('%d (%.0f%%)', trxs_vs_peaks_median, 100*trxs_vs_peaks_median/total_peaks)

# 842 (12%)
trxs_vs_bkg_median <- as.integer(median(trxs_vs_bkg_gg_df$n_signif))
sprintf('%d (%.0f%%)', trxs_vs_bkg_median, 100*trxs_vs_bkg_median/total_peaks)

###

rand_vs_peaks_gg_df <- filter(peak_gg_df, startsWith(RNA, 'rand'))
rand_vs_bkg_gg_df   <- filter(bkg_gg_df,  startsWith(RNA, 'rand'))

# 2677 (39%)
rand_vs_peaks_median <- as.integer(median(rand_vs_peaks_gg_df$n_signif))
sprintf('%d (%.0f%%)', rand_vs_peaks_median, 100*rand_vs_peaks_median/total_peaks)

# 645 (9%)
rand_vs_bkg_median <- as.integer(median(rand_vs_bkg_gg_df$n_signif))
sprintf('%d (%.0f%%)', rand_vs_bkg_median, 100*rand_vs_bkg_median/total_peaks)

