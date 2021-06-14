library(here)
library(vegan)
source(here("03_import_seq_data.R"))

# subset anemones
anemones <- subset_samples(phy, SampleType == "Anemone")

# generate rarefaction curves
rarecurve(t(otu_table(anemones)),
          step = 600,
          sample = 12000,
          ylab = "Sample ASVs",
          xlab = "Sample reads",
          main = "Gnotobiotic samples")

rm(list = ls())
