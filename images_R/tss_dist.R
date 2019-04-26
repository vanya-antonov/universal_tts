
source("lib.R")

###

get_dist_group_df <- function(dist_mtx, min_dist, max_dist)
{
  # min_dist, max_dist are in kb
  group_name <- sprintf("%.0f - %.0f kb", min_dist/1000, max_dist/1000)
  data.frame(
    dist_group = rep(group_name, ncol(dist_mtx)),
    dist_freq = apply(min_dist <= dist_mtx & dist_mtx < max_dist, 2, sum) / nrow(dist_mtx))
}

get_rand_gg_df <- function(rand_dist_fn)
{
  bkg_dist_mtx <- paste0(DATA_DIR, rand_dist_fn) %>%
    read.delim(as.is=TRUE, header=FALSE) %>%
    as.matrix() %>%
    abs()
  
  rbind(
    get_dist_group_df(bkg_dist_mtx, 0,    1000),
    get_dist_group_df(bkg_dist_mtx, 1000, 2000),
    get_dist_group_df(bkg_dist_mtx, 2000, 3000),
    get_dist_group_df(bkg_dist_mtx, 3000, 4000),
    get_dist_group_df(bkg_dist_mtx, 4000, 5000))
}

get_uni_gg_df <- function(uni_dist_fn)
{
  uni_dist_v <- paste0(DATA_DIR, uni_dist_fn) %>%
    read.delim(as.is=TRUE, header=FALSE) %>%
    pull(11) %>%
    abs()
  
  # Create 1-col mtx
  uni_dist_mtx <- matrix(uni_dist_v, nrow = length(uni_dist_v), ncol = 1)
  rbind(
    get_dist_group_df(uni_dist_mtx, 0,    1000),
    get_dist_group_df(uni_dist_mtx, 1000, 2000),
    get_dist_group_df(uni_dist_mtx, 2000, 3000),
    get_dist_group_df(uni_dist_mtx, 3000, 4000),
    get_dist_group_df(uni_dist_mtx, 4000, 5000))
  
}

get_pvalue_df <- function(uni_gg_df, rand_gg_df)
{
  # Create list of vectors: https://stackoverflow.com/a/38542871/310453
  rand_lv <- split(rand_gg_df$dist_freq, uni_gg_df$dist_group)
  uni_lv <- split(uni_gg_df$dist_freq, uni_gg_df$dist_group)
  
  empirial_pvalue <- function(name){ sum(rand_lv[[name]] >= uni_lv[[name]]) / length(rand_lv[[name]]) }
  pvalue_v <- sapply(uni_gg_df$dist_group, empirial_pvalue)
  
  data.frame(dist_group = uni_gg_df$dist_group, pvalue = pvalue_v) %>%
    mutate(stars = case_when(pvalue <= 0.01 ~ "**",
                             pvalue <= 0.05 ~ "*",
                             pvalue <= 1 ~ ""))
}

get_tss_dist_ggplot <- function(base_name, title='')
{
  uni_dist_fn <- sprintf('tss_dist.%s.txt', base_name)
  rand_dist_fn <- sprintf('tss_dist.%s.rand.txt', base_name)
  
  uni_gg_df <- get_uni_gg_df(uni_dist_fn)
  rand_gg_df <- get_rand_gg_df(rand_dist_fn)
  pvalue_df <- get_pvalue_df(uni_gg_df, rand_gg_df)

  num_uni <- paste0(DATA_DIR, uni_dist_fn) %>%
    read.delim(as.is=TRUE, header=FALSE) %>%
    nrow()
  subT <- paste0("N = ", num_uni)

  ggplot(rand_gg_df) +
    aes(x = factor(dist_group), y = dist_freq) +
    geom_boxplot(outlier.shape = NA, fill = 'gray') +
    geom_point(data = uni_gg_df, shape=21, size = 4, fill = 'red') +
    geom_text(data = pvalue_df, aes(label = stars), y = 0.42, size = 10) +
    # Add the '%' sign to the Y-axis: https://stackoverflow.com/a/41098629/310453
    scale_y_continuous(
      labels = scales::percent_format(),
      limits = c(0, 0.45)) +
    xlab('Absolute distance to TSS') +
    ylab('% of genomic regions') +
    ggtitle(title, subtitle = subT)
}

chop_gg <- get_tss_dist_ggplot('uni_tts_chop_seq', 'Universal ChOP-seq peaks')
bkg_gg <- get_tss_dist_ggplot('uni_tts_background', 'Universal Background regions')
capture_gg <- get_tss_dist_ggplot('uni_tts_capture_seq', 'Universal Capture-seq peaks')

plot_grid(
  chop_gg, bkg_gg, capture_gg,
  ncol = 1,
  labels = 'AUTO'
)
ggsave('tss_dist.pdf', path = OUT_DIR, width = 8, height = 11)
