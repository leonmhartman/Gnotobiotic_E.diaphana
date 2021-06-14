library(here)
library(decontam)
library(microbiome)

source(here("03_import_seq_data.R"))

# identify contaminants
consList <- isContaminant(seqtab = phy, neg = "neg", method = "prevalence")

# get contaminant names
cons <- rownames(consList)[consList$contaminant=="TRUE"]


# to get info on the contaminants, uncomment the following code
# and run it ON THE FILE WITH THE CONTAMINANT ASVs IN-SITU, then
# combine the consPer.csv and taxonomy.csv file data

# # subset the anemone samples FROM THE PRE-DECONTAM PHYLOSEQ FILE to find out
# # what percentage the contaminants were in the original anemone samples
# vvv <- subset_samples(phy, SampleType == "Anemone")
# # merge the samples
# yyy <- merge_samples(vvv, group = "SampleType", fun = sum)
# # transform counts to percentages
# yyy <- transform_sample_counts(yyy, function(x) 100 * x/sum(x))
# # extract the cons percentage data
# zzz <- prune_taxa(x = yyy, taxa = cons)
# # write otu table to dataframe
# xxx <- data.frame(t(zzz@otu_table))
# # write xxx to csv
# write.csv(here(x = xxx, row.names = TRUE, file = "contamAnemones.csv"))
# # subset the contaminant ASVs
# gnotoPhyCons <- prune_taxa(phy, taxa = cons)
# # write the contaminants to a file for reference (when required i.e. not EVERY time this script is run)
# contaminants <- write_phyloseq(here(x = gnotoPhyCons, type = 'TAXONOMY'))
#
# # subset the artemia samples FROM THE PRE-DECONTAM PHYLOSEQ FILE to find out
# # what percentage the contaminants were in the original artemia samples
# vvv <- subset_samples(phy, SampleType == "Artemia")
# # merge the samples
# yyy <- merge_samples(vvv, group = "SampleType", fun = sum)
# # transform counts to percentages
# yyy <- transform_sample_counts(yyy, function(x) 100 * x/sum(x))
# # extract the cons percentage data
# zzz <- prune_taxa(x = yyy, taxa = cons)
# # write otu table to dataframe
# xxx <- data.frame(t(zzz@otu_table))
# # write xxx to csv
# write.csv(here(x = xxx, row.names = TRUE, file = "contamArtemia.csv"))
#
# rm(vvv)
# rm(xxx)
# rm(yyy)
# rm(zzz)
# rm(gnotoPhyCons)
# rm(contaminants)


#                                   Anemone Artemia
# 32980c3819b9e36c51a512dd11770683	0.0025  0
# 6d9b5102954749bab307469fd1bc861c	0.0007  0
# a43cc15d35c5f4cf5db945c8b532c680	0.0041  0
# aecb3e24bb57cb642935054dc6f8fbd3	0.0081  0.0051
# e7dcbd80845ef8b0a075b50c17cc150c	0.0003  0
# 79fea7073a32159fbf592c02ee449582	0.0601  0
# de90abf55c878c685dfa474a85660dfe	0.0033  0

# 32980c3819b9e36c51a512dd11770683	Proteobacteria	Gammaproteobacteria				
# 6d9b5102954749bab307469fd1bc861c	Proteobacteria	Alphaproteobacteria	Rhizobiales	          Rhizobiaceae	      Mesorhizobium	
# a43cc15d35c5f4cf5db945c8b532c680	Proteobacteria	Alphaproteobacteria	Rhizobiales	          Xanthobacteraceae	  Bradyrhizobium	
# aecb3e24bb57cb642935054dc6f8fbd3	Proteobacteria	Alphaproteobacteria	Rhizobiales	          Xanthobacteraceae	  Afipia	          Bradyrhizobium sp. CCH10-C7
# e7dcbd80845ef8b0a075b50c17cc150c	Bacteroidetes	  Bacteroidia	        Sphingobacteriales	  Sphingobacteriaceae	Sphingobacterium	
# 79fea7073a32159fbf592c02ee449582	Proteobacteria	Gammaproteobacteria	Betaproteobacteriales	Burkholderiaceae	  Ralstonia	
# de90abf55c878c685dfa474a85660dfe  Proteobacteria	Gammaproteobacteria	Betaproteobacteriales	Burkholderiaceae		


# remove the contaminants from the main phyloseq file
phy <- remove_taxa(phy, taxa = cons)
# 4620 ASVs across all samples

rm(cons)
rm(consList)
