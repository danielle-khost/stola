library(tidyverse)

syn <- read.table("~/github/toga4cnee/tests/test_runs/pigeonChickenStrict/overlap_stats.synteny.flank1000.tsv", sep="\t", header=T)

syn <- syn %>%
  mutate(flank_seq_identity = as.numeric(replace(flank_seq_identity, flank_seq_identity == ".", "0")))

ggplot(syn, aes(x=flank_seq_identity, fill=orthology_status)) +
  geom_histogram() +
  facet_wrap(vars(orthology_status), ncol=1)

ggplot(syn, aes(x=flank_fraction_gap, fill=orthology_status)) +
  geom_histogram() +
  facet_wrap(vars(orthology_status), ncol=1)

ggplot(syn, aes(x=flank_seq_identity,y=flank_fraction_gap))+
  geom_bin2d(bins=10) +
  scale_fill_continuous(type = "viridis") +
  theme_classic() +
  facet_wrap(vars(orthology_status), ncol=1)
