library(tidyverse)

syn <- read.table("~/github/toga4cnee/tests/test_runs/colLivi/overlap_stats.synteny.tsv", sep="\t", header=T)

 syn <- syn %>%
   mutate(distance_to_nearest_cnee_neighbor = as.numeric(replace(distance_to_nearest_cnee_neighbor, distance_to_nearest_cnee_neighbor == ".", NA)))

#2D histogram
ggplot(syn, aes(x=flank_seq_identity,y=flank_fraction_gap))+
  geom_bin2d(bins=10) +
  scale_fill_continuous(type = "viridis") +
  theme_bw() +
  facet_grid(orthology_status ~ flank_size)

#density 
synSub <- syn %>% subset(flank_size == 100) 

##distance to cnee
ggplot(synSub, aes(x=distance_to_nearest_cnee_neighbor,color=orthology_status)) +
  geom_density() +
  theme_bw() +
  scale_x_log10() +
  facet_grid(vars(orthology_status))

syn %>% subset(flank_size == 100 & orthology_status == "lost") %>%
  ggplot(., aes(x=distance_to_nearest_aligned_block)) +
  geom_density() +
  theme_bw() +
  scale_x_log10()
