
source("lib.R")

###

get_coverage_gg <- function(fn)
{
  #fn <- 'coverage_chop_uni.txt'
  #fn <- 'coverage_chop_other.txt'
  # fn <- 'coverage_other_bkg.txt'
  df <- read.delim(paste0(DATA_DIR, fn), as.is=TRUE, header=TRUE)
  title <- sprintf("%s:%s (%d bp)", df[1, 'name'], df[1, 'strand'], nrow(df))
  df %>%
    mutate(frac_rna = num_rna / 153,
           nt_type = ifelse(poly_ga == 1, 'poly_ga', 'other_nt')) %>%
           #nt_type = factor(nt_type, levels=c('other_nt', 'poly_ga'))) %>%
    ggplot() +
    aes(x = coord) +
#    geom_area(aes(y = frac_rna)) +
    geom_area(aes(y = num_rna)) +
    geom_rug(aes(col = nt_type), sides = 'b') +
    # scale_y_continuous(
    #   labels = scales::percent,   # https://stackoverflow.com/a/41098629/310453
    #   limits = c(0, 1)) +   
    ylim(0, 153) +
    ggtitle(title) +
    xlab("DNA region position (bp)") +
    # ylab("% Query RNAs that form triplexes") +
    ylab("Number of interacting RNAs") +
    scale_colour_manual(
      name = "",
      values = c('white', 'red'),
      breaks = c('other_nt', 'poly_ga'),
      labels = c('', "Poly-GA element")) +
    theme(legend.position="bottom")  # Move legend to the bottom
}

W = 9
H = 6
get_coverage_gg('coverage_chop_uni.txt')
ggsave('coverage_chop_uni.pdf', path = OUT_DIR, width = W, height = H)

get_coverage_gg('coverage_bkg_uni.txt')
ggsave('coverage_bkg_uni.pdf', path = OUT_DIR, width = W, height = H)

get_coverage_gg('coverage_capture_uni.txt')
ggsave('coverage_capture_uni.pdf', path = OUT_DIR, width = W, height = H)

