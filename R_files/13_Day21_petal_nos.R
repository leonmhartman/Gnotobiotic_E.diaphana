library(here)

source(here("08_outlier_deletion.R"))

# get the numbers of common and unique ASVs in the Day 21 control and treated anemones

# subset Day 21 samples
C1 <- subset_samples(anemones, sampNames == "gc211")
C2 <- subset_samples(anemones, sampNames == "gc212")
C3 <- subset_samples(anemones, sampNames == "gc213")
C4 <- subset_samples(anemones, sampNames == "gc214")
C5 <- subset_samples(anemones, sampNames == "gc215")
C6 <- subset_samples(anemones, sampNames == "gc216")
T1 <- subset_samples(anemones, sampNames == "gt211")
T2 <- subset_samples(anemones, sampNames == "gt212")
T3 <- subset_samples(anemones, sampNames == "gt213")
T4 <- subset_samples(anemones, sampNames == "gt214")
T6 <- subset_samples(anemones, sampNames == "gt216")

# remove zero sum ASVs
C1 <- prune_taxa((taxa_sums(C1) > 0), C1)
C2 <- prune_taxa((taxa_sums(C2) > 0), C2)
C3 <- prune_taxa((taxa_sums(C3) > 0), C3)
C4 <- prune_taxa((taxa_sums(C4) > 0), C4)
C5 <- prune_taxa((taxa_sums(C5) > 0), C5)
C6 <- prune_taxa((taxa_sums(C6) > 0), C6)
T1 <- prune_taxa((taxa_sums(T1) > 0), T1)
T2 <- prune_taxa((taxa_sums(T2) > 0), T2)
T3 <- prune_taxa((taxa_sums(T3) > 0), T3)
T4 <- prune_taxa((taxa_sums(T4) > 0), T4)
T6 <- prune_taxa((taxa_sums(T6) > 0), T6)

# create 'not in' operator
'%!in%' <-  Negate('%in%')

# get number of ASVs unique to each sample
sum(taxa_names(C1) %!in% unique(c(taxa_names(C2),taxa_names(C3),taxa_names(C4),taxa_names(C5),taxa_names(C6))))
sum(taxa_names(C2) %!in% unique(c(taxa_names(C1),taxa_names(C3),taxa_names(C4),taxa_names(C5),taxa_names(C6))))
sum(taxa_names(C3) %!in% unique(c(taxa_names(C2),taxa_names(C1),taxa_names(C4),taxa_names(C5),taxa_names(C6))))
sum(taxa_names(C4) %!in% unique(c(taxa_names(C2),taxa_names(C3),taxa_names(C1),taxa_names(C5),taxa_names(C6))))
sum(taxa_names(C5) %!in% unique(c(taxa_names(C2),taxa_names(C3),taxa_names(C4),taxa_names(C1),taxa_names(C6))))
sum(taxa_names(C6) %!in% unique(c(taxa_names(C2),taxa_names(C3),taxa_names(C4),taxa_names(C5),taxa_names(C1))))
sum(taxa_names(T1) %!in% unique(c(taxa_names(T2),taxa_names(T3),taxa_names(T4),taxa_names(T6))))
sum(taxa_names(T2) %!in% unique(c(taxa_names(T1),taxa_names(T3),taxa_names(T4),taxa_names(T6))))
sum(taxa_names(T3) %!in% unique(c(taxa_names(T2),taxa_names(T1),taxa_names(T4),taxa_names(T6))))
sum(taxa_names(T4) %!in% unique(c(taxa_names(T2),taxa_names(T3),taxa_names(T1),taxa_names(T6))))
sum(taxa_names(T6) %!in% unique(c(taxa_names(T2),taxa_names(T3),taxa_names(T4),taxa_names(T1))))

# get number of ASVs common to all samples
length(Reduce(intersect, list(taxa_names(C1),taxa_names(C2),taxa_names(C3),taxa_names(C4),taxa_names(C5),taxa_names(C6))))
length(Reduce(intersect, list(taxa_names(T1),taxa_names(T2),taxa_names(T3),taxa_names(T4),taxa_names(T6))))

# get total number of ASVs in each group
phyC21 <- subset_samples(anemones, SampleType == "Anemone" & Treatment == "Control" & Day == "21")
phyC21 <- prune_taxa((taxa_sums(phyC21) > 0), phyC21)
nrow(phyC21@tax_table)
phyT21 <- subset_samples(anemones, SampleType == "Anemone" & Treatment == "Treatment" & Day == "21")
phyT21 <- prune_taxa((taxa_sums(phyT21) > 0), phyT21)
nrow(phyT21@tax_table)

rm(list = ls())
