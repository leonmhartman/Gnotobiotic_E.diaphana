library(here)
library(ggplot2)

source(here("08_outlier_deletion.R"))

# encode days as factor
sample_data(anemones)$Day <- factor(sample_data(anemones)$Day)

# day 0 and 21 anemones
day0day21 <- subset_samples(anemones, Day == "0" | Day == "21")

# prune taxa
day0day21 <- prune_taxa((taxa_sums(day0day21) > 0), day0day21)

day0day21.nmds <- ordinate(day0day21, method = "NMDS", distance = "bray")

plot_ordination(
  day0day21, 
  ordination = day0day21.nmds, 
  type = "samples",
  color = "Treatment", 
  shape = "Day"
  ) + 
  theme_bw() + 
  geom_point(size = 5) +
  stat_ellipse()

rm(list = ls())
