library(phyloseq)

# import feature table
maireTable <- read.table("maire_table.tsv",
                         header = TRUE,
                         sep = "\t",
                         row.names = 1)

# import taxonomy table
maireTax <- read.table("maire_tax.tsv",
                       sep = "\t",
                       header = FALSE,
                       fill = TRUE,
                       row.names = 1)

# add taxonomy levels
colnames(maireTax) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# import metadata
maireMeta <- read.table('maire_meta.txt',
                        header = TRUE,
                        sep = '\t',
                        row.names = 1)

# fix column name with hyphen
colnames(maireMeta)[1] <- "sampleType"

# convert data frames to matrices
maireTable_mat <- as.matrix(maireTable)
maireTax_mat <- as.matrix(maireTax)

# combine feature, taxonomy, and metadata tables into a phyloseq object
mairePhy <- phyloseq(otu_table(maireTable_mat, taxa_are_rows = T),
                tax_table(maireTax_mat),
                sample_data(maireMeta))

# remove zero-sum ASVs
mairePhy <- prune_taxa((taxa_sums(mairePhy) > 0), mairePhy)
# 484 ASVs

rm(maireTable_mat)
rm(maireTax_mat)
