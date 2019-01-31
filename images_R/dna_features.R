
library(dplyr)
library(tidyr)  # gather()
library(ggplot2)
library(cowplot)   # To have ggplots side-by-side: plot_grid()

###

DATA_DIR <- "../data/"
OUT_DIR <- "../images/"

theme_set(theme_bw(base_size = 19))  # increase the font size: https://stackoverflow.com/a/11955412/310453

###

read_all_data <- function(gc_fn, pur_fn, uni_tts_fn)
{
  # gc_fn <- 'background.len_gc'
  # pur_fn <- 'background.len_purine'
  # uni_tts_fn <- 'uni_tts_background_hg38.bed'
  gc_df <- read.delim(paste0(DATA_DIR, gc_fn), as.is=TRUE, header=FALSE,
                      col.names = c('dna', 'len', 'gc'))
  pur_df <- read.delim(paste0(DATA_DIR, pur_fn), as.is=TRUE, header=FALSE,
                       col.names = c('dna', 'len', 'purine', 'poly_purine'))
  uni_tts_v <- read.delim(paste0(DATA_DIR, uni_tts_fn), as.is=TRUE, header=FALSE)$V4
  inner_join(gc_df, pur_df, by = c('dna', 'len')) %>%
    mutate(uni_tts = dna %in% uni_tts_v)
}

chop_df <- read_all_data('chop_seq.len_gc', 'chop_seq.len_purine', 'uni_tts_chop_seq_hg38.bed')
bkg_df <- read_all_data('background.len_gc', 'background.len_purine', 'uni_tts_chop_seq_hg38.bed')

chop_df <- mutate(chop_df, type = ifelse(uni_tts, 'uni_tts', 'other'))
bkg_df$type <- 'bkg'

feat_names <- c('gc' = 'GC content',
                'purine' = 'Purine content',
                'poly_purine' = 'Poly-purine content')

dna_names <- c('bkg' = sprintf('Control DNA regions (n = %d)', nrow(bkg_df)),
               'uni_tts' = sprintf('Universal ChOP-seq peaks (n = %d)', sum(chop_df$uni_tts)),
               'other' = sprintf('Other ChOP-seq peaks (n = %d)', sum(!chop_df$uni_tts)))

bind_rows(chop_df, bkg_df) %>%
  gather('gc', 'purine', 'poly_purine', key = 'feature', value = 'value') %>%
  # Define the order of the sub-plots (facets)
  mutate(feature = factor(feature, levels = names(feat_names))) %>%
  # Define the order of the boxplots
  mutate(type = factor(type, levels = names(dna_names))) %>%
  ggplot() +
  aes(x = type, y = value/100, fill = type) +
  geom_boxplot(outlier.shape = NA) +
  # geom_jitter(width = 0.2, alpha = 0.5, size = 0.5) +
  facet_grid(. ~ feature, labeller = as_labeller(feat_names)) +
  theme(legend.position="top") +
  # Use custom colors and provide the legend description without modifying the underlying factor
  # http://www.cookbook-r.com/Graphs/Legends_(ggplot2)/#with-fill-and-color
  scale_fill_manual(name=NULL,  # Remove legend title
                    values = c('white', 'red', 'gray'),
                    breaks=names(dna_names),
                    labels=dna_names) +
  # Legend items should be in one column: https://stackoverflow.com/a/18400725/310453
  guides(fill=guide_legend(ncol=1)) +
  # Add the '%' sign to the Y-axis: https://stackoverflow.com/a/41098629/310453
  scale_y_continuous(labels = scales::percent_format()) +
  # Remove X-labels: https://stackoverflow.com/a/35090981/310453
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())
ggsave('dna_features.pdf', path=OUT_DIR, width = 9, height = 6)

