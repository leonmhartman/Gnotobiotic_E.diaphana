library(here)
library(dplyr)
library(plotrix)

source(here("08_outlier_deletion.R"))

# normalise by rarefying to compensate for differences in sample number and sequencing effort
anemones <- rarefy_even_depth(anemones, sample.size = min(sort(sample_sums(anemones))), rngseed = 1)

# subset
anemones <- subset_samples(anemones, SampleType != "Artemia")
anemones <- prune_taxa((taxa_sums(anemones) > 0), anemones)

# add grouping variable
variable1 = as.character(get_variable(anemones, "Treatment"))
variable2 = as.character(get_variable(anemones, "Day"))
sample_data(anemones)$groups <- mapply(paste0, variable1, variable2, collapse = "_")

# merge control and treated counts by day
merged <- merge_samples(anemones, "groups")

# transform absolute counts to relative abundance
percents <- transform_sample_counts(merged, function(x) 100 * x/sum(x))

# identify the 20 most abundant ASVs
top20 <- which(colnames(percents@otu_table) %in% top_taxa(percents, n=20))

# merge other ASVs
others <- 1:ncol(percents@otu_table)
others <- others[-top20]
archOthers <- colnames(percents@otu_table)[(others[1])]
top20 <- merge_taxa(percents, others, archetype = archOthers)

# melt data
q <- psmelt(top20)

# change stacking so least abundant are at top
q$OTU <- reorder(q$OTU, -q$Abundance)
q$OTU <- factor(q$OTU, levels = rev(levels(q$OTU)))

# set x-axis order
q$Sample <- factor(q$Sample, levels = c("Control0","Control1","Control3",
                                        "Control7","Control14","Control21",
                                        "Treatment0","Treatment1","Treatment3",
                                        "Treatment7","Treatment14","Treatment21"))

# specify bar colours
barCols <- c("#2317FF", "#D01B00", "#00A80E", "#F48BFF", "#FFD617", "#00E7B9", "#D045FF",
             "#8B3300", "#FF4817", "#BEAFF5", "#706F49", "#EF7C6E", "#444FC7", "#A3FF2E",
             "#A06D46", "#5A0007", "#DBBF69", "#59EDA5", "#0B6000", "#7B45FF", "#070C05",
             "#2317FF", "#D01B00", "#00A80E", "#F48BFF", "#FFD617", "#00E7B9", "#D045FF",
             "#8B3300", "#FF4817", "#BEAFF5", "#706F49", "#EF7C6E", "#444FC7", "#A3FF2E")

# generate ggplot2 object
ggplot(q, aes_string(x = "Sample", y = "Abundance", fill = "OTU")) +  
  geom_bar(stat = 'identity', position = "stack") +
  ylab("Reads assigned to ASV (%)") +
  theme_bw() +
  scale_fill_manual(values = barCols) +
  guides(fill = guide_legend(ncol = 1)) +
  theme(axis.text.x = element_text(angle = 90))

# write data to file
write.table(q, file='bar_chart_legend_data.tsv', sep='\t')

rm(list = ls())