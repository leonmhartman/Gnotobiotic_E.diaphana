library(here)
source(here("06_import_seq_data_bh_adj.R"))

# remove two outlier samples
anemones <- subset_samples(anemones, sampNames != "gt03" & sampNames != "gt11")

# prune taxa
anemones <- prune_taxa((taxa_sums(anemones) > 0), anemones)
