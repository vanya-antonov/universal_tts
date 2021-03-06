
source("lib.R")   # get_log_pvalue_mtx(), find_uni_tts()

###

peak_init_df <- read.table(paste0(DATA_DIR, "pvalue_chop_seq.txt.gz"), header=TRUE, row.names=1)
bkg_init_df <- read.table(paste0(DATA_DIR, "pvalue_background.txt.gz"), header=TRUE, row.names=1)
capture_init_df <- read.table(paste0(DATA_DIR, "pvalue_capture_seq.txt.gz"), header=TRUE, row.names=1)

# Are all the p-values valid? (must be TRUE)
all(0 <= peak_init_df & peak_init_df <= 1)
all(0 <= bkg_init_df & bkg_init_df <= 1)
all(0 <= capture_init_df & capture_init_df <= 1)

peak_data <- get_log_pvalue_mtx(peak_init_df)
bkg_data <- get_log_pvalue_mtx(bkg_init_df)
capture_data <- get_log_pvalue_mtx(capture_init_df)

peak_gg_df <- find_uni_tts(peak_data) 
bkg_gg_df <- find_uni_tts(bkg_data) 
capture_gg_df <- find_uni_tts(capture_data) 

get_q_fraction_gg <- function(gg_df)
{
  gg_df <- gg_df[order(-gg_df$q_fraction),]
  gg_df$x <- 1:nrow(gg_df)
  type_descr <- c(
    sprintf('Universal TTS \n(n = %d)', sum(gg_df$type == '1_uni_tts')),
    sprintf('Other regions \n(n = %d)', sum(gg_df$type == '2_other')))
  ggplot(gg_df) +
    aes(x = x, y = q_fraction, fill = type) +
    geom_area() +
    scale_fill_manual(labels = type_descr, values = c("red", "black")) +
    xlab('DNA region index') +
    ylab(paste0('% query RNAs with\nadj p-value < ', PVALUE_THR)) +
    theme(legend.position = "top") +
    theme(legend.key.size = unit(2, 'lines')) +
    theme(legend.title=element_blank())
}

plot_uni_tts_chop <- get_q_fraction_gg(peak_gg_df) + ggtitle('ChOP-seq peaks')
plot_uni_tts_bkg <- get_q_fraction_gg(bkg_gg_df) + ggtitle('Control DNA regions')
plot_uni_tts_capture <- get_q_fraction_gg(capture_gg_df) + ggtitle('Shared Capture-seq peaks')


# PLOT ALL! ----
plot_grid(
  plot_uni_tts_chop, plot_uni_tts_bkg, plot_uni_tts_capture,
  ncol = 1, labels = 'AUTO')
ggsave(paste0('uni_tts.pdf'), path = OUT_DIR, width = 8, height = 14)


# Export all uni_tts ----
save_uni_tts_as_bed <- function(gg_df, out_fn)
{
  # gg_df <- peak_gg_df
  # out_fn <- "uni_tts_chop_hg38.bed"
  filter(gg_df, type == "1_uni_tts") %>%
    separate(dna, into = c('chr', 'start', 'end'),
             sep = '[:-]', remove = FALSE) %>%
    select(chr, start, end, dna) %>%
    write_tsv(paste0(DATA_DIR, out_fn), col_names = FALSE)
}
save_uni_tts_as_bed(peak_gg_df, "uni_tts_chop_seq_hg38.bed")
save_uni_tts_as_bed(bkg_gg_df, "uni_tts_background_hg38.bed")
save_uni_tts_as_bed(capture_gg_df, "uni_tts_capture_seq_hg38.bed")


# Statistics ----
total_peaks_chop <- nrow(peak_gg_df)
total_peaks_capture <- nrow(capture_gg_df)

uni_tts_peak_n <- filter(peak_gg_df, type == "1_uni_tts") %>% nrow()
uni_tts_bkg_n <- filter(bkg_gg_df, type == "1_uni_tts") %>% nrow()
uni_tts_capture_n <- filter(capture_gg_df, type == "1_uni_tts") %>% nrow()

# 1107 (16.3%)
sprintf('%d (%.1f%%)', uni_tts_peak_n, 100*uni_tts_peak_n/total_peaks_chop)

# 41 (0.6%)
sprintf('%d (%.1f%%)', uni_tts_bkg_n, 100*uni_tts_bkg_n/total_peaks_chop)

# 2151 (40.0%)
sprintf('%d (%.1f%%)', uni_tts_capture_n, 100*uni_tts_capture_n/total_peaks_capture)
