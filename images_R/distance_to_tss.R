
source("lib.R")

###

DIST_COLS <- c('chr', 'start', 'end', 'name', 'dist')

read_dist_data <- function(dist_fn, uni_tts_fn, dna_group)
{
  # dist_fn <- 'chop_seq_tss_dist.txt'
  # uni_tts_fn <- 'uni_tts_chop_seq_hg38.bed'
  # dna_group <- 'chop'
  uni_tts_v <- read.delim(paste0(DATA_DIR, uni_tts_fn), as.is=TRUE, header=FALSE)$V4
  paste0(DATA_DIR, dist_fn) %>%
    read.delim(as.is=TRUE, header=FALSE, col.names = DIST_COLS) %>%
    # Create groups for historgam
    # mutate(dist_group = case_when(abs(dist) <= 3000 ~ "< 3 kb",
    #                               abs(dist) > 3000 ~ "> 3 kb")) %>%
    mutate(dist_group = case_when(abs(dist) <= 1000 ~ "0 - 1 kb",
                                  abs(dist) <= 2000 ~ "1 - 2 kb",
                                  abs(dist) <= 3000 ~ "2 - 3 kb",
                                  abs(dist) <= 4000 ~ "3 - 4 kb",
                                  abs(dist) <= 5000 ~ "4 - 5 kb",
                                  abs(dist) > 5000 ~ ">5 kb")) %>%
    # Add info about uni_tts
    mutate(dna_type = ifelse(name %in% uni_tts_v, 'uni', 'other')) %>%
    mutate(dna_group = dna_group,
           group_type = paste0(dna_group, '_', dna_type))
    
}

chop_df <- read_dist_data('chop_seq_tss_dist.txt', 'uni_tts_chop_seq_hg38.bed', 'chop')
bkg_df <- read_dist_data('background_tss_dist.txt', 'uni_tts_background_hg38.bed', 'bkg')
capture_df <- read_dist_data('capture_seq_tss_dist.txt', 'uni_tts_capture_seq_hg38.bed', 'capture')

total_v <- c(
  'chop' = paste0(DATA_DIR, 'chop_seq_hg38.bed') %>% read.delim(header=FALSE) %>% nrow(),
  'bkg' = paste0(DATA_DIR, 'background_hg38.bed') %>% read.delim(header=FALSE) %>% nrow(),
  'capture' = paste0(DATA_DIR, 'capture_seq_hg38.bed') %>% read.delim(header=FALSE) %>% nrow(),
  'chop_uni' = paste0(DATA_DIR, 'uni_tts_chop_seq_hg38.bed') %>% read.delim(header=FALSE) %>% nrow(),
  'bkg_uni' = paste0(DATA_DIR, 'uni_tts_background_hg38.bed') %>% read.delim(header=FALSE) %>% nrow(),
  'capture_uni' = paste0(DATA_DIR, 'uni_tts_capture_seq_hg38.bed') %>% read.delim(header=FALSE) %>% nrow())
total_v['chop_other'] = total_v['chop'] - total_v['chop_uni']
total_v['bkg_other'] = total_v['bkg'] - total_v['bkg_uni']
total_v['capture_other'] = total_v['capture'] - total_v['capture_uni']

stats_df <- bind_rows(bkg_df, chop_df, capture_df) %>%
  group_by(dna_group, dna_type, group_type, dist_group) %>%
  summarise(num_dna = n_distinct(name)) %>%
  mutate(dist_freq = num_dna / total_v[group_type]) %>%
  ungroup()

stats_df %>% select(dist_group) %>% distinct()
type_names = c('other' = 'Other regions',
               'uni' = 'Universal TTSs')

group_names <- c( 'chop' = 'ChOP-seq',
                  'bkg' = 'Control',
                  'capture' = 'Capture-seq')

dist_names <- c(
  '0 - 1 kb', '1 - 2 kb', '2 - 3 kb', '3 - 4 kb', '4 - 5 kb', '>5 kb'
)

stats_df %>%
  # Define the order of DNA groups (columns)
  mutate(dna_group = factor(dna_group, levels = names(group_names))) %>%
  # Define the order of boxplots (DNA types)
  mutate(dna_type = factor(dna_type, levels = names(type_names))) %>%
  # Define the order of Distances
  mutate(dist_group = factor(dist_group, levels = dist_names)) %>%
  ggplot() +
  aes(x = dist_group, y = dist_freq, fill = dna_type) +
  geom_bar(col = 'black', stat = 'identity', position = 'dodge', width = 0.8) +
  # Add the '%' sign to the Y-axis: https://stackoverflow.com/a/41098629/310453
  scale_y_continuous(labels = scales::percent_format()) +
  facet_grid(. ~ dna_group, labeller = as_labeller(group_names)) +
  # Move legend to the top
  theme(legend.position="top") +
  # http://www.cookbook-r.com/Graphs/Legends_(ggplot2)/#with-fill-and-color
  scale_fill_manual(name=NULL,                       # Remove legend title
                    values = c('gray', 'red'),      # Use custom colors
                    breaks=rev(names(type_names)),   # Reverse the items in legend
                    labels=type_names) +             # provide the legend description
  xlab('Absolute distance to TSS') +
  ylab('% of genomic regions') +
  # Rotate X-axis labels 45 degrees
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave('distance_to_tss.pdf', path = OUT_DIR, width = 11, height = 7)
