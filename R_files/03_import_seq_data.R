library(here)
library(phyloseq)

# import feature table
gnotoTable <- read.table(here("gnoto_table.tsv"),
                         header = TRUE,
                         sep = "\t",
                         row.names = 1)

# import taxonomy table
gnotoTax <- read.table(here("gnoto_tax.tsv"),
                       sep = "\t",
                       header = FALSE,
                       fill = TRUE,
                       row.names = 1)

# add taxonomy levels
colnames(gnotoTax) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# import metadata
gnotoMeta <- read.table(here('gnoto_meta.txt'),
                        header = TRUE,
                        sep = '\t',
                        row.names = 1)

# fix sample name with hyphen
row.names(gnotoMeta)[76] <- "ge.blank"

# add values for decontam
gnotoMeta$neg <- FALSE
gnotoMeta$neg[76:79] <- TRUE

# convert data frames to matrices
gnotoTable_mat <- as.matrix(gnotoTable)
gnotoTax_mat <- as.matrix(gnotoTax)

# combine feature, taxonomy, and metadata tables into a phyloseq object
phy <- phyloseq(otu_table(gnotoTable_mat, taxa_are_rows = T),
                tax_table(gnotoTax_mat),
                sample_data(gnotoMeta))

# # this code adds the QIIME-generated tree to the phyloseq object, but since the
# # tree wasn't used for any of the analyses in this work, I didn't add it.
# gnotoTree <- read_tree(here("gnoto_tree.nwk"))
# phy <- merge_phyloseq(phy, gnotoTree)

# remove zero-sum ASVs
phy <- prune_taxa((taxa_sums(phy) > 0), phy)
# 4627 ASVs

rm(gnotoTable_mat)
rm(gnotoTax_mat)
