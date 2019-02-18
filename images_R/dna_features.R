
source("lib.R")

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
    mutate(dna_type = ifelse(dna %in% uni_tts_v, 'uni_tts', 'other'))
  # mutate(uni_tts = dna %in% uni_tts_v)
}

chop_df <- read_all_data('chop_seq.len_gc', 'chop_seq.len_purine', 'uni_tts_chop_seq_hg38.bed')
chop_df$dna_group <- 'chop'

bkg_df <- read_all_data('background.len_gc', 'background.len_purine', 'uni_tts_background_hg38.bed')
bkg_df$dna_group <- 'bkg'

capture_df <- read_all_data('capture_seq.len_gc', 'capture_seq.len_purine', 'uni_tts_capture_seq_hg38.bed')
capture_df$dna_group <- 'capture'

type_names = c('other' = 'Other regions',
               'uni_tts' = 'Universal TTSs')

feat_names <- c('gc' = 'GC content',
                'purine' = 'GA content',
                'poly_purine' = 'Poly-GA')

group_names <- c( 'chop' = 'ChOP-seq',
                  'bkg' = 'Control',
                  'capture' = 'Capture-seq')

all_labels <- as_labeller(c(group_names, feat_names))

bind_rows(chop_df, bkg_df, capture_df) %>%
  mutate(poly_purine_len = len * poly_purine / 100) %>%
  #gather('len', 'gc', 'purine', 'poly_purine', 'poly_purine_len', key = 'feature', value = 'value') %>%
  gather('gc', 'purine', 'poly_purine', key = 'feature', value = 'value') %>%
  # Define the order of features (rows)
  mutate(feature = factor(feature, levels = names(feat_names))) %>%
  # Define the order of DNA groups (columns)
  mutate(dna_group = factor(dna_group, levels = names(group_names))) %>%
  # Define the order of boxplots (DNA types)
  mutate(dna_type = factor(dna_type, levels = names(type_names))) %>%
  ggplot() +
  aes(x = dna_type, y = value/100, fill = dna_type) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(feature ~ dna_group, labeller = all_labels) +
  # Add the '%' sign to the Y-axis: https://stackoverflow.com/a/41098629/310453
  scale_y_continuous(labels = scales::percent_format()) +
  # Move legend to the top
  theme(legend.position="top") +
  # http://www.cookbook-r.com/Graphs/Legends_(ggplot2)/#with-fill-and-color
  scale_fill_manual(name=NULL,                       # Remove legend title
                    values = c('white', 'red'),      # Use custom colors
                    breaks=rev(names(type_names)),   # Reverse the items in legend
                    labels=type_names) +             # provide the legend description
  # Remove axis labels: https://stackoverflow.com/a/35090981/310453
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())
ggsave('dna_features.pdf', path=OUT_DIR, width = 8, height = 8)


# statistics ----
# 1 bkg                3.43
# 2 other              8.7 
# 3 uni_tts           35.6 
bind_rows(chop_df, bkg_df) %>%
  group_by(dna_type) %>%
  summarise(median_poly_pur = median(poly_purine))

