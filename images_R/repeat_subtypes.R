
source("lib.R")

###

REP_COLS <- c('chr', 'start', 'end', 'name', 'rep_chr', 'rep_start', 'rep_end', 'type_subtype', 'ovlp_len')

read_repeat_data <- function(repeat_fn, uni_tts_fn, dna_group)
{
  # repeat_fn <- 'chop_seq_repeats.txt'
  # uni_tts_fn <- 'uni_tts_chop_seq_hg38.bed'
  # dna_group <- 'chop'
  uni_tts_v <- read.delim(paste0(DATA_DIR, uni_tts_fn), as.is=TRUE, header=FALSE)$V4
  paste0(DATA_DIR, repeat_fn) %>%
    read.delim(as.is=TRUE, header=FALSE, col.names = REP_COLS) %>%
    separate(type_subtype, sep = "::", into = c('rep_type', 'rep_subtype'), remove = FALSE) %>%
    filter(rep_type %in% c("Simple_repeat", "Low_complexity")) %>%
    group_by(name, type_subtype, rep_type, rep_subtype) %>%
    summarise(ovlp_sum_len = sum(ovlp_len)) %>%
    # There must be at least 10 bp of the overlap inside the region
    filter(ovlp_sum_len > 10) %>%
    # Add info about uni_tts
    mutate(dna_type = ifelse(name %in% uni_tts_v, 'uni', 'other')) %>%
    mutate(dna_group = dna_group,
           group_type = paste0(dna_group, '_', dna_type))
}

chop_df <- read_repeat_data('chop_seq_repeats.txt', 'uni_tts_chop_seq_hg38.bed', 'chop')
bkg_df <- read_repeat_data('background_repeats.txt', 'uni_tts_background_hg38.bed', 'bkg')
capture_df <- read_repeat_data('capture_seq_repeats.txt', 'uni_tts_capture_seq_hg38.bed', 'capture')

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
  group_by(dna_group, dna_type, group_type, type_subtype, rep_type, rep_subtype) %>%
  summarise(dna_with_rep = n_distinct(name)) %>%
  mutate(rep_freq = dna_with_rep / total_v[group_type]) %>%
  ungroup()

frequent_subtypes <- stats_df %>%
  filter(rep_freq >= 0.05) %>%
  arrange(rep_type, -rep_freq) %>%
  pull(type_subtype) %>%
  unique()

group_type_names <- c(
  'chop_uni' = 'Universal ChOP peaks ',
  'bkg_uni' = 'Universal Control regions ',
  'capture_uni' = 'Universal Capture peaks ',
  'chop_other' = 'Other ChOP peaks',
  'bkg_other' = 'Other Control regions',
  'capture_other' = 'Other Capture peaks')

stats_df %>%
  # Select frequent subtypes only
  filter(type_subtype %in% frequent_subtypes) %>%
  # Add rows with zeros for missin repeat subtypes in some groups
  spread(type_subtype, rep_freq, fill = 0) %>%
  gather(frequent_subtypes, key='type_subtype', value='rep_freq') %>%
  # Define the order of bars in each repeat (group_type)
  mutate(group_type = factor(group_type, levels = group_type_names %>% names() %>% rev())) %>%
  # Define the order of repeat subtypes
  mutate(type_subtype = factor(type_subtype, levels = rev(frequent_subtypes))) %>%
  ggplot() +
  aes(x = type_subtype, y = rep_freq, fill = group_type) +
  geom_bar(col = 'black', stat = 'identity', position = 'dodge', width = 0.8) +
  # Add the '%' sign to the Y-axis: https://stackoverflow.com/a/41098629/310453
  scale_y_continuous(labels = scales::percent_format()) +
  coord_flip() +
  theme(axis.title.y=element_blank(),
        legend.title=element_blank(),
        legend.position="top") +
  # Custom colors: http://www.cookbook-r.com/Graphs/Legends_(ggplot2)/#with-fill-and-color
  scale_fill_manual(values = rev(c('red', 'orange', 'yellow', 'blue4', 'blue', 'deepskyblue')),
                    #values = rev(c('red', 'orange', 'gold', 'black', 'gray40', 'gray80')),
                    breaks=rev(names(group_type_names)),   # Reverse the items in legend
                    labels=group_type_names) +             # provide the legend description
  #xlab("Repeat type") +
  ylab("% of genomic regions") +
  # About reverse = TRUE: https://stackoverflow.com/a/49448745/310453
  # About ncol = 2: https://stackoverflow.com/a/36087339/310453
  guides(fill = guide_legend(reverse = TRUE, ncol = 2))
ggsave('repeat_subtypes.pdf', path = OUT_DIR, width = 9, height = 11)


