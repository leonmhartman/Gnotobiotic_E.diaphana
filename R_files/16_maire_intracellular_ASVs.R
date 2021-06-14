library(here)
library(phyloseq)

source(here("15_maire_import_seq_data.R"))

# import tolerant ASVs
tolerantASVs <- read.table("tolerant_ASVs.tsv",
                         header = TRUE,
                         sep = "\t",
                         row.names = 1)

# subset "intracellular" samples
intra <- subset_samples(mairePhy, sampleType == "intra")
intra <- prune_taxa((taxa_sums(intra) > 0), intra)
# 265 intracellular ASVs

# subset "loosely-associated" samples
loose <- subset_samples(mairePhy, sampleType == "loose")
loose <- prune_taxa((taxa_sums(loose) > 0), loose)
# 168 intracellular ASVs

# subset "closely-associated" samples
close <- subset_samples(mairePhy, sampleType == "close")
close <- prune_taxa((taxa_sums(close) > 0), close)
# 282 intracellular ASVs

# which tolerant ASVs are associated with each Maire sample type?
intersect(rownames(tolerantASVs), rownames(intra@tax_table))
# [1] "613a6844484c3b2f90d38399d2698624" "e3c9fb8d8e882e9f6d98ff1d7c32f03f"
# [3] "ae8a6381f4e78d50e3fad0e2faad8ea0" "e0cc6a95596c1068ec99eb67cea8d93e"
# [5] "53818a706e38c1584f139d2f90fbd8df" "acd191a5f307d0579b3126166f102435"
# [7] "26a930ff8d7ae487ae9b3beaf07b53ca"
intersect(rownames(tolerantASVs), rownames(loose@tax_table))
# [1] "ae8a6381f4e78d50e3fad0e2faad8ea0" "e0cc6a95596c1068ec99eb67cea8d93e"
# [3] "acd191a5f307d0579b3126166f102435"
intersect(rownames(tolerantASVs), rownames(close@tax_table))
# [1] "613a6844484c3b2f90d38399d2698624" "e3c9fb8d8e882e9f6d98ff1d7c32f03f"
# [3] "ae8a6381f4e78d50e3fad0e2faad8ea0" "e0cc6a95596c1068ec99eb67cea8d93e"
# [5] "53818a706e38c1584f139d2f90fbd8df" "acd191a5f307d0579b3126166f102435"
# [7] "26a930ff8d7ae487ae9b3beaf07b53ca"

# 7 tolerant ASVs were found in the Maire intracellular samples.
# The same 7 ASVs were also found in the loosely-associated samples.
# 3 tolerant ASVs were also found in the closely-associated samples.

rm(list = ls())
