source(here("06_import_seq_data_bh_adj.R"))

# encode days as factor
sample_data(anemones)$Day <- factor(sample_data(anemones)$Day)

# create plot file
anemones.nmds <- ordinate(anemones, method = "NMDS", distance = "bray")

plot_ordination(
  anemones,
  ordination = anemones.nmds, 
  type = "samples",
  color = "Treatment", 
  shape = "Day",
  label = "sampNames"
) + 
  theme_bw() + 
  geom_point(size = 5) +
  stat_ellipse()

rm(list = ls())
