
library(dplyr)
library(ggplot2)
library(cowplot)  # To have ggplots side-by-side: plot_grid()

###

OUT_DIR <- "../images/"

theme_set(theme_bw(base_size = 15))  # increase the font size: https://stackoverflow.com/a/11955412/310453

###

# Number of triplexes ----
DATA_DIR <- "~/Projects/2018/Yulia/DataLog/1119.F1000.rand_logtpot/data/"

all_data <- list.files(DATA_DIR, pattern=".*\\.tpx$", full.names=TRUE) %>%
  lapply(read.table, as.is=TRUE, header = TRUE) %>%
  Reduce(rbind, .)

title <- sprintf("MEG3 length = %s nt", unique(all_data$RNA_len))
subT <- "Parameters: -l 10 -e 10 -fr off"
plot_num_tpx_vs_dna_len <- ggplot(all_data) +
  aes(x = factor(DNA_len), y = Num_tpx) +
  geom_boxplot(outlier.size = 0.5) +
  xlab("Random DNA length") +
  ylab("Number of triplexes") +
  ggtitle(title, subtitle = subT)


# lambda ----
DATA_DIR <- "~/2019/2019/Yulia/DataLog/0110.F1000.poisson/data/"

fn <- paste0(DATA_DIR, '_RAND_100_vs_100.tpx.gz')
init_data <- read.table(fn, as.is=TRUE, header = TRUE)
head(init_data)

all_lambda <- init_data %>%
  group_by(RNA_len, DNA_len) %>%
  summarise(lambda = mean(Num_tpx_abs))

lambda_lm <- lm(lambda ~ RNA_len + DNA_len, data = all_lambda)
summary(lambda_lm)
coef(lambda_lm)
str(summary(lambda_lm))

fixed_RNA_len <- 1500
plot_lambda_lm <- all_lambda %>%
  filter(RNA_len == fixed_RNA_len) %>%
  ggplot() +
  aes(x = DNA_len, y = lambda) +
  geom_point() +
  geom_abline(intercept = coef(lambda_lm)[1] + fixed_RNA_len * coef(lambda_lm)[2],
              slope = coef(lambda_lm)[3],
              linetype = "dashed") +
  ylim(0, NA) +
  ggtitle(paste0("Random RNA length = ", fixed_RNA_len),
          subtitle = sprintf("Adjusted R-squared = %.3f", summary(lambda_lm)$adj.r.squared)) +
  xlab("Random DNA length (bp)") +
  ylab("Average number\nof triplexes")


# Histograms with Poisson ----
fixed_DNA_len <- c(200, 400, 800, 1600)
gg_pois <- expand.grid(Num_tpx_abs = seq(0, 15), DNA_len = fixed_DNA_len) %>%
  mutate(RNA_len = fixed_RNA_len,
         lambda = predict(lambda_lm, data.frame(RNA_len, DNA_len)),
         pois_value = dpois(Num_tpx_abs, lambda))

plot_num_tpx_with_pred_pois <- init_data %>%
  filter(RNA_len == fixed_RNA_len, DNA_len %in% fixed_DNA_len) %>%
  ggplot() +
  aes(x = Num_tpx_abs) +
  # Use freqs in geom_bar: https://stackoverflow.com/a/49372165/310453
  geom_bar(aes(y = ..prop..)) +
  geom_path(aes(y = pois_value), data = gg_pois, colour="red") +
  geom_point(aes(y = pois_value), data = gg_pois, colour="red") +
  facet_grid(. ~ DNA_len) +
  xlim(NA, 15) +
  ggtitle(paste0("Random RNA length = ", fixed_RNA_len)) +
  xlab("Number of predicted triplexes") +
  ylab("Frequency")


# PLOT ALL! ----
row1 <- plot_grid(
  plot_num_tpx_vs_dna_len, plot_lambda_lm,
  nrow = 1, labels=c('a', 'b'))

row2 <- plot_grid(
  plot_num_tpx_with_pred_pois,
  nrow = 1, labels=c('c'))

plot_grid(
  row1, row2,
  rel_heights = c(1, 0.9),
  ncol = 1)
ggsave('rand_dna.pdf', path=OUT_DIR, width = 9, height = 6)

