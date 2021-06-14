#!/bin/sh

#---------------------------------------------------------#

# Import raw sequencing files (refer to maire_mainfest.txt for list of
# files to download from NCBI BioProject ID PRJNA650221)

rm -r ~/maire/seqs/demux
mkdir ~/maire/seqs/demux

qiime tools import \
--type SampleData[PairedEndSequencesWithQuality] \
--input-path ~/maire/maire_qiime2_manifest.txt \
--input-format PairedEndFastqManifestPhred33V2 \
--output-path ~/maire/seqs/demux/demuxed.qza

#---------------------------------------------------------#

rm -r ~/maire/seqs/trim

qiime cutadapt trim-paired \
--i-demultiplexed-sequences ~/maire/seqs/demux/demuxed.qza \
--p-front-f AGGATTAGATACCCTGGTA \
--p-front-r CRRCACGAGCTGACGAC \
--p-error-rate 0.20 \
--p-minimum-length 100 \
--p-cores 5 \
--output-dir ~/maire/seqs/trim/ \
--verbose

#---------------------------------------------------------#

# Generate viewable file to check quality

rm ~/maire/seqs/trim/trimmed_sequences.qzv

qiime demux summarize \
--i-data ~/maire/seqs/trim/trimmed_sequences.qza \
--o-visualization ~/maire/seqs/trim/trimmed_sequences.qzv \
--verbose

# View on https://view.qiime2.org

#---------------------------------------------------------#

# Perform dada2 error correction, and trimming

rm -r ~/maire/seqs/dada2out

qiime dada2 denoise-paired \
--i-demultiplexed-seqs ~/maire/seqs/trim/trimmed_sequences.qza \
--p-trunc-len-f 196 \
--p-trunc-len-r 120 \
--p-n-threads 0 \
--output-dir ~/maire/seqs/dada2out \
--verbose


# Rename output files

mv ~/maire/seqs/dada2out/table.qza \
~/maire/seqs/dada2out/maire_table.qza
mv ~/maire/seqs/dada2out/representative_sequences.qza \
~/maire/seqs/dada2out/maire_rep_seqs.qza
mv ~/maire/seqs/dada2out/denoising_stats.qza \
~/maire/seqs/dada2out/maire_denoising_stats.qza

#---------------------------------------------------------#

# Generate viewable files

rm ~/maire/seqs/dada2out/maire_table.qzv
rm ~/maire/seqs/dada2out/maire_rep_seqs.qzv
rm ~/maire/seqs/dada2out/maire_denoising_stats.qzv

qiime feature-table summarize \
--i-table ~/maire/seqs/dada2out/maire_table.qza \
--m-sample-metadata-file ~/maire/maire_meta.txt \
--o-visualization ~/maire/seqs/dada2out/maire_table.qzv \
--verbose

qiime feature-table tabulate-seqs \
--i-data ~/maire/seqs/dada2out/maire_rep_seqs.qza \
--o-visualization ~/maire/seqs/dada2out/maire_rep_seqs.qzv \
--verbose

qiime metadata tabulate \
--m-input-file ~/maire/seqs/dada2out/maire_denoising_stats.qza \
--o-visualization ~/maire/seqs/dada2out/maire_denoising_stats.qzv \
--verbose

# View on https://view.qiime2.org

#---------------------------------------------------------#

# Assign taxonomy with SILVA v132 database

rm -r ~/maire/taxonomy

qiime feature-classifier classify-sklearn \
--i-classifier ~/maire/silva_132_16s_v5v6_classifier.qza \
--i-reads ~/maire/seqs/dada2out/maire_rep_seqs.qza \
--p-n-jobs -2 \
--output-dir ~/maire/taxonomy/ \
--verbose


# Rename output file

mv ~/maire/taxonomy/classification.qza \
~/maire/taxonomy/maire_taxonomy.qza

#---------------------------------------------------------#

# Generate viewable summary file of taxonomic assignments

rm ~/maire/taxonomy/maire_taxonomy.qzv

qiime metadata tabulate \
--m-input-file ~/maire/taxonomy/maire_taxonomy.qza \
--o-visualization ~/maire/taxonomy/maire_taxonomy.qzv \
--verbose

# View on https://view.qiime2.org

#---------------------------------------------------------#

# Remove unassigned, mitochondria and chloroplast features, if they exist

rm ~/maire/seqs/maire_table_filtered.qza

qiime taxa filter-table \
--i-table ~/maire/seqs/dada2out/maire_table.qza \
--i-taxonomy ~/maire/taxonomy/maire_taxonomy.qza \
--p-exclude Unassigned,Mitochondria,Chloroplast \
--o-filtered-table ~/maire/seqs/dada2out/maire_table_filtered.qza \
--verbose

#---------------------------------------------------------#

# Generate viewable file of filtered table

rm ~/maire/seqs/dada2out/maire_table.qzv

qiime feature-table summarize \
--i-table ~/maire/seqs/dada2out/maire_table_filtered.qza \
--m-sample-metadata-file ~/maire/maire_meta.txt \
--o-visualization ~/maire/seqs/dada2out/maire_table_filtered.qzv \
--verbose

# View on https://view.qiime2.org

#---------------------------------------------------------#

# Perform alignment on rep seqs

rm -r ~/maire/tree
mkdir ~/maire/tree

qiime alignment mafft \
--i-sequences ~/maire/seqs/dada2out/maire_rep_seqs.qza \
--p-n-threads 5 \
--o-alignment ~/maire/tree/aligned_maire_rep_seqs.qza \
--verbose

#--------------------------------------------------#

# Mask highly variable regions of the alignment

rm ~/maire/tree/masked_aligned_maire_rep_seqs.qza

qiime alignment mask \
--i-alignment ~/maire/tree/aligned_maire_rep_seqs.qza \
--o-masked-alignment ~/maire/tree/masked_aligned_maire_rep_seqs.qza \
--verbose

#--------------------------------------------------#

# Generate phylogenetic tree

rm ~/maire/tree/maire_unrooted_tree.qza

qiime phylogeny fasttree \
--i-alignment ~/maire/tree/masked_aligned_maire_rep_seqs.qza \
--p-n-threads 1 \
--o-tree ~/maire/tree/maire_unrooted_tree.qza \
--verbose

#--------------------------------------------------#

# Apply mid-point rooting to tree

rm ~/maire/tree/maire_rooted_tree.qza

qiime phylogeny midpoint-root \
--i-tree ~/maire/tree/maire_unrooted_tree.qza \
--o-rooted-tree ~/maire/tree/maire_rooted_tree.qza \
--verbose

#--------------------------------------------------#

# Output files for R

mkdir ~/maire/output

qiime tools export \
--input-path ~/maire/tree/maire_unrooted_tree.qza \
--output-path ~/maire/output/

qiime tools export \
--input-path ~/maire/taxonomy/maire_taxonomy.qza \
--output-path ~/maire/output/

qiime tools export \
--input-path ~/maire/seqs/dada2out/maire_table_filtered.qza \
--output-path ~/maire/output/

biom convert -i ~/maire/output/feature-table.biom \
--to-tsv \
-o ~/maire/output/maire_table.tsv

#--------------------------------------------------#

# Edit output files

# Rename tree file
mv ~/maire/output/tree.nwk ~/maire/output/maire_tree.nwk

# Remove biom file
rm ~/maire/output/feature-table.biom

# Remove header from feature table # Ubuntu coding is: sed -i '1d'
sed -i '' '1d' ~/maire/output/maire_table.tsv

# Remove header from taxonomy table # Ubuntu coding is: sed -i '1,2d'
sed -i '' '1,2d' ~/maire/output/taxonomy.tsv

# Remove confidence values from taxonomy table
cut -f 1,2 ~/maire/output/taxonomy.tsv > ~/maire/output/maire_tax.tsv

# Remove old taxonomy file
rm ~/maire/output/taxonomy.tsv

# Replace semi-colons in taxonomy table with tabs # Ubuntu coding is: sed -i 's/;/\t/g'
# For OSX, "<CTRL-V><TAB>" must be entered at the control line
sed -i '' 's/;/<CTRL-V><TAB>/g' ~/maire/output/maire_tax.tsv

#--------------------------------------------------#

# The end
