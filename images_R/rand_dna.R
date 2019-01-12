
library(MASS)     # fitdistr()
library(dplyr)
library(ggplot2)
library(cowplot)  # To have ggplots side-by-side: plot_grid()

###

OUT_DIR <- "../images/"
#DATA_DIR <- "~/Projects_2018/2018/Yulia/DataLog/1211.F1000.poisson/data/"

theme_set(theme_bw(base_size = 15))  # increase the font size: https://stackoverflow.com/a/11955412/310453

###

# Number of triplexes ----
DATA_DIR <- "~/Projects_2018/2018/Yulia/DataLog/1119.F1000.rand_logtpot/data/"

all_data <- list.files(DATA_DIR, pattern=".*\\.tpx$", full.names=TRUE) %>%
  lapply(read.table, as.is=TRUE, header = TRUE) %>%
  Reduce(rbind, .)

title <- sprintf("RNA length = %s", unique(all_data$RNA_len))
subT <- "Parameters: -l 10 -e 10 -fr off"
plot_num_tpx_vs_dna_len <- ggplot(all_data) +
  aes(x = factor(DNA_len), y = Num_tpx) +
  geom_boxplot(outlier.size = 0.5) +
  xlab("DNA length") +
  ylab("Number of triplexes") +
  ggtitle(title, subtitle = subT)


# MEG3_peak_len_hist ----
DATA_DIR <- "~/Projects_2018/2018/Yulia/DataLog/0518.F1000.presentation/data/"
meg_peaks <- read.delim(paste0(DATA_DIR, '_MEG3.peaks_wo_N.len_purine'),
                        as.is=TRUE, header=FALSE, col.names = c('name','len','purine'))

title <- sprintf("MEG3 peak lengths (n=%d)", nrow(meg_peaks))
subT <- sprintf("median length = %.1f bp", median(meg_peaks$len))
plot_peak_len_hist <- ggplot(meg_peaks) +
  aes(x = len) +
  geom_histogram(binwidth = 30, col = "black", fill = "white") +
  xlim(0, 1500) +
  labs(x = "Peak length (bp)", y = "Number of peaks") +
  ggtitle(title, subtitle = subT)


# lambda ----
DATA_DIR <- "~/Projects_2018/2018/Yulia/DataLog/1211.F1000.poisson/data/"

init_data <- read.table(paste0(DATA_DIR, '_MEG3_vs_RAND_DNA.tpx.gz'), as.is=TRUE, header = TRUE)
#head(init_data)

all_lambda <- init_data %>%
  group_by(DNA_len) %>%
  summarise(lambda = fitdistr(Num_tpx_abs, densfun="poisson")$estimate)

lambda_lm <- lm(lambda ~ DNA_len, data = all_lambda)

pred_pois <- expand.grid(Num_tpx_abs = seq(0, 20), lambda = all_lambda$lambda) %>%
  right_join(all_lambda, by = "lambda") %>%
  mutate(pred_lambda = predict(lambda_lm, data.frame(DNA_len = DNA_len)),
         pred_pois_value = dpois(Num_tpx_abs, lambda = pred_lambda))


plot_lambda_lm <- ggplot(all_lambda) +
  aes(x = DNA_len, y = lambda) +
  geom_point() +
  geom_abline(intercept = coef(lambda_lm)[1],
              slope = coef(lambda_lm)[2],
              linetype = "dashed") +
  ggtitle("MEG3 vs random DNA",
          subtitle=sprintf("lambda = %s + %s*DNA_len",
                  format(coef(lambda_lm)[1], digits = 2),
                  format(coef(lambda_lm)[2], digits = 2))) +
  xlab("Random DNA length (bp)") +
  ylab("Average number of\npredicted triplexes")

plot_num_tpx_with_pred_pois <- ggplot(init_data) +
  aes(x = Num_tpx_abs) +
  # Use freqs in geom_bar: https://stackoverflow.com/a/49372165/310453
  geom_bar(aes(y = ..prop..)) +
  geom_path(aes(y = pred_pois_value), data = pred_pois, colour="red") +
#  geom_point(aes(y = pred_pois_value), data = pred_pois, colour="red") +
  xlab("Number of predicted triplexes") +
  ylab("Frequency") +
  facet_grid(. ~ DNA_len)


# PLOT ALL! ----
row1 <- plot_grid(
  plot_num_tpx_vs_dna_len, NULL,
  nrow = 1, labels=c('a', ''))

row2 <- plot_grid(
  plot_peak_len_hist, plot_lambda_lm,
  #rel_widths = c(1, 1.3),
  nrow = 1, labels=c('b', 'c'))

row3 <- plot_grid(
  plot_num_tpx_with_pred_pois,
  nrow = 1, labels=c('d'))

plot_grid(
  row1, row2, row3,
  rel_heights = c(1, 1, 0.9),
  ncol = 1)
ggsave('rand_dna.pdf', path=OUT_DIR, width = 9, height = 9)

