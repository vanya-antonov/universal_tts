
library(dplyr)
library(ggplot2)

###

OUT_DIR <- "../images/"
DATA_DIR <- "../data/"

theme_set(theme_bw(base_size = 19))  # increase the font size: https://stackoverflow.com/a/11955412/310453

###

info <- read.delim(paste0(DATA_DIR, 'expressed_trxs.txt'), as.is=TRUE, header=FALSE)
colnames(info) <- c('trx', 'name', 'len', 'gc', 'reads', 'rpkm')

info <- subset(info, rpkm > 1)
very_similar <- subset(info, 1400 < len & len < 1800 & 55 < gc & gc < 60)
other_trxs <- info[!rownames(info) %in% rownames(very_similar), ]

# MEG3: NR_002766	1595	57.55
meg3_df <- data.frame(gc=57.55, len=1595)

ggplot(other_trxs) +
  aes(x = gc, y = len) +
  geom_point(alpha = 0.2, size = 1) +
  geom_point(data = very_similar, col = "red", alpha = 0.5, size = 1) +
#  geom_point(data = meg3_df, col = "yellow", size = 5) +
  geom_point(data = meg3_df, col = "green", size = 5) +
  ylim(0, 4000) +
  xlim(20, 80) +
  xlab("Transcript GC content (%)") +
  ylab("Transcript length (nt)") +
  ggtitle(sprintf("Total number of transcripts with RPKM > %.0f: %d", min(info$rpkm), nrow(info)),
          subtitle = sprintf("Number of selected transcripts: %d", nrow(very_similar)))

# Increase dpi for better resolution: https://stackoverflow.com/q/47222764/310453
ggsave("expressed_trxs.png", path=OUT_DIR, dpi = 150)
