library(here)
library(ggplot2)
library(mixOmics)

source(here("08_outlier_deletion.R"))

# encode days as factor
sample_data(anemones)$Day <- factor(sample_data(anemones)$Day)

# - - - - - - - - - - - - - - - - - - #

# nMDS - all days

# create plot file
alldays.nmds <- ordinate(anemones, method = "NMDS", distance = "bray")

plot_ordination(
  anemones, 
  ordination = alldays.nmds, 
  type = "samples",
  color = "Treatment", 
  shape = "Day"
) + 
  theme_bw() + 
  geom_point(size = 5) +
  stat_ellipse() +
  scale_x_reverse()

# export plot
ggsave("nMDS_all_days.pdf", scale = 2, width = 9, height = 6, units = "cm")

# - - - - - - - - - - - - - - - - - - #

# nMDS - day 0 vs day 21

# subset data
day0day21 <- subset_samples(anemones, Day == "0" | Day == "21")

# prune taxa
day0day21 <- prune_taxa((taxa_sums(day0day21) > 0), day0day21)

# create plot file
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
  stat_ellipse() +
  scale_shape_manual(values=c(16, 8)) +
  scale_x_reverse()
  
# export plot
ggsave("nMDS_day0day21.pdf", scale = 2, width = 9, height = 6, units = "cm")

# - - - - - - - - - - - - - - - - - - #

# PCA - all days

# get ASV data
asv <- as.data.frame(t(anemones@otu_table))

# add small non-zero value to allow CLR transformation
asv <- asv + 0.001

# The PCA is generated using the mixOmics package and methodology: http://mixomics.org/
# Under this methodology, low count data is often filtered. However, because low count
# ASVs were the primary drivers of change in this experiment, they have been retained.

# get metadata
meta <- as.data.frame(anemones@sam_data)

# add group i.d.
meta$Group <- paste0(meta$Treatment, "_", meta$Day)

# apply transformation
data.clr <- logratio.transfo(asv, logratio = 'CLR')

# create PCA file
pca.clr <- pca(data.clr)

# pch to match nMDS
pch_lookup <- c("Control_0" = 16,
                "Control_1" = 17,
                "Control_14"= 7,
                "Control_21" = 8,
                "Control_3" = 15,
                "Control_7" = 3,
                "Treatment_0" = 16,
                "Treatment_1" = 17,
                "Treatment_14" = 7,
                "Treatment_21" = 8,
                "Treatment_3" = 15,
                "Treatment_7" = 3)

# plot PCA
plotIndiv(pca.clr, group = meta$Group, size.title = 0, ellipse = TRUE, pch = pch_lookup, legend = FALSE)

# export plot
ggsave("PCA_all_days.pdf", scale = 2, width = 7, height = 6, units = "cm")

# - - - - - - - - - - - - - - - - - - #

# PCA - day 0 vs day 21

# subset data
day0day21 <- subset_samples(anemones, Day == "0" | Day == "21")

# prune taxa
day0day21 <- prune_taxa((taxa_sums(day0day21) > 0), day0day21)

# get ASV data
asv <- as.data.frame(t(day0day21@otu_table))

# add small non-zero value to allow CLR transformation
asv <- asv + 0.001

# The PCA is generated using the mixOmics package and methodology: http://mixomics.org/
# Under this methodology, low count data is often filtered. However, because low count
# ASVs were the primary drivers of change in this experiment, they have been retained.

# get metadata
meta <- as.data.frame(day0day21@sam_data)

# add group i.d.
meta$Group <- paste0(meta$Treatment, "_", meta$Day)

# apply transformation
data.clr <- logratio.transfo(asv, logratio = 'CLR')

# create PCA file
pca.clr <- pca(data.clr)

# pch to match nMDS
pch_lookup <- c("Control_0" = 16,
                "Control_21" = 8,
                "Treatment_0" = 16,
                "Treatment_21" = 8)

# plot PCA
plotIndiv(pca.clr, group = meta$Group, size.title = 0, ellipse = TRUE, pch = pch_lookup, legend = FALSE)

# export plot
ggsave("PCA_day0day21.pdf", scale = 2, width = 7, height = 6, units = "cm")

rm(list = ls())
