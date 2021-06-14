library(here)
library(phyloseq)

source(here("05_decontam.R"))

# subset anemones
anemones <- subset_samples(phy, SampleType == "Anemone")

# prune taxa
anemones <- prune_taxa((taxa_sums(anemones) > 0), anemones)

# import B/H data
bh <- read.table(here("ddPCR_counts.txt"),
                 header = TRUE,
                 sep = "\t")

# calculate BH ratios (values were corrected by subtracting NTC background signal)
bh$ratio <- bh$Bcorrected/bh$Hcorrected

# get BH ratios
BH <- bh[2:73,c(1,9)]

# change order of BH ratios to match order of otu table
BH <- BH[order(match(BH$sample, colnames(anemones@otu_table))),]

# scale counts according to BH ratios
otu_table(anemones) <- t(round(t(otu_table(anemones)) * BH$ratio * 1000))
# scaling is by multiplication using the BH ratio, as per Jian 2020.
# this produces '16S reads per host cell x 10^3' which is easy to interpret and accounts for
# sample size differences.
# note that the data is also multiplied by 1000 because the BH ratio values are 0-1 so many
# counts end up being <1, which creates issues with some downstream analyses, such as the
# alpha diversity comparison which includes a rarefying step (without the extra multiplication,
# the rarefying value is so low, many ASVs are discarded).

# remove gt215

# add column of sample names for labelling
anemones@sam_data$sampNames <- rownames(anemones@sam_data)

# remove gt215
anemones <- subset_samples(anemones, sampNames != "gt215")

# prune taxa
anemones <- prune_taxa((taxa_sums(anemones) > 0), anemones)

rm(bh)
rm(BH)
