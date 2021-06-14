#!/bin/sh

# Execution time = 55 mins on 2018 Mac Mini, 6-core i7, 32 Gb RAM

#---------------------------------------------------------#

# Import raw sequencing files

mkdir ~/gnoto/seqs/demux

qiime tools import \
--type SampleData[PairedEndSequencesWithQuality] \
--input-path ~/gnoto/seqs/raw \
--input-format CasavaOneEightSingleLanePerSampleDirFmt \
--output-path ~/gnoto/seqs/demux/gnoto_demuxed.qza

#---------------------------------------------------------#

# Generate viewable file to determine truncation lengths

qiime demux summarize \
--i-data ~/gnoto/seqs/demux/gnoto_demuxed.qza \
--o-visualization ~/gnoto/seqs/demux/gnoto_demuxed.qzv

# View on https://view.qiime2.org

#---------------------------------------------------------#

# Perform dada2 error correction and trimming

qiime dada2 denoise-paired \
--i-demultiplexed-seqs ~/gnoto/seqs/demux/gnoto_demuxed.qza \
--p-trunc-len-f 250 \
--p-trunc-len-r 199 \
--p-trim-left-f 19 \
--p-trim-left-r 17 \
--p-n-threads 0 \
--output-dir ~/gnoto/seqs/dada2out \
--verbose


# Rename output files

mv ~/gnoto/seqs/dada2out/table.qza \
~/gnoto/seqs/dada2out/gnoto_table.qza
mv ~/gnoto/seqs/dada2out/representative_sequences.qza \
~/gnoto/seqs/dada2out/gnoto_rep_seqs.qza
mv ~/gnoto/seqs/dada2out/denoising_stats.qza \
~/gnoto/seqs/dada2out/gnoto_denoising_stats.qza

#---------------------------------------------------------#

# Generate viewable files

qiime feature-table summarize \
--i-table ~/gnoto/seqs/dada2out/gnoto_table.qza \
--m-sample-metadata-file ~/gnoto/gnoto_qiime2_map.txt \
--o-visualization ~/gnoto/seqs/dada2out/gnoto_table.qzv \
--verbose

qiime feature-table tabulate-seqs \
--i-data ~/gnoto/seqs/dada2out/gnoto_rep_seqs.qza \
--o-visualization ~/gnoto/seqs/dada2out/gnoto_rep_seqs.qzv \
--verbose

qiime metadata tabulate \
--m-input-file ~/gnoto/seqs/dada2out/gnoto_denoising_stats.qza \
--o-visualization ~/gnoto/seqs/dada2out/gnoto_denoising_stats.qzv \
--verbose

# View on https://view.qiime2.org

#---------------------------------------------------------#

# Remove features >=290 nt and generate viewable files

qiime feature-table filter-seqs \
--i-data ~/gnoto/seqs/dada2out/gnoto_rep_seqs.qza \
--m-metadata-file ~/gnoto/seqs/dada2out/gnoto_rep_seqs.qza \
--p-where 'length(sequence) < 290' \
--o-filtered-data ~/gnoto/seqs/dada2out/gnoto_rep_seqs_filtered.qza

qiime feature-table filter-features \
--i-table ~/gnoto/seqs/dada2out/gnoto_table.qza \
--m-metadata-file ~/gnoto/seqs/dada2out/gnoto_rep_seqs.qza \
--p-where 'length(sequence) < 290' \
--o-filtered-table ~/gnoto/seqs/dada2out/gnoto_table_filtered.qza


qiime feature-table tabulate-seqs \
--i-data ~/gnoto/seqs/dada2out/gnoto_rep_seqs_filtered.qza \
--o-visualization ~/gnoto/seqs/dada2out/gnoto_rep_seqs_filtered.qzv \
--verbose

qiime feature-table summarize \
--i-table ~/gnoto/seqs/dada2out/gnoto_table_filtered.qza \
--m-sample-metadata-file ~/gnoto/gnoto_qiime2_map.txt \
--o-visualization ~/gnoto/seqs/dada2out/gnoto_table_filtered.qzv \
--verbose

# View on https://view.qiime2.org

#---------------------------------------------------------#

# Assign taxonomy with SILVA v132 database trained with
# 784 + 1061 primers

qiime feature-classifier classify-sklearn \
--i-classifier ~/gnoto/silva_132_16s_v5v6_classifier.qza \
--i-reads ~/gnoto/seqs/dada2out/gnoto_rep_seqs_filtered.qza \
--p-n-jobs -2 \
--output-dir ~/gnoto/taxonomy/ \
--verbose


# Rename output file

mv ~/gnoto/taxonomy/classification.qza \
~/gnoto/taxonomy/gnoto_taxonomy.qza

#---------------------------------------------------------#

# Generate viewable summary file of taxonomic assignments

qiime metadata tabulate \
--m-input-file ~/gnoto/taxonomy/gnoto_taxonomy.qza \
--o-visualization ~/gnoto/taxonomy/gnoto_taxonomy.qzv \
--verbose

# View on https://view.qiime2.org

#---------------------------------------------------------#

# Remove unassigned, mitochondria and chloroplast features, if they exist

qiime taxa filter-table \
--i-table ~/gnoto/seqs/dada2out/gnoto_table_filtered.qza \
--i-taxonomy ~/gnoto/taxonomy/gnoto_taxonomy.qza \
--p-exclude Unassigned,Mitochondria,Chloroplast \
--o-filtered-table ~/gnoto/seqs/dada2out/gnoto_table_filtered.qza \
--verbose

#---------------------------------------------------------#

# Generate viewable file of filtered feature table

qiime feature-table summarize \
--i-table ~/gnoto/seqs/dada2out/gnoto_table_filtered.qza \
--m-sample-metadata-file ~/gnoto/gnoto_qiime2_map.txt \
--o-visualization ~/gnoto/seqs/dada2out/gnoto_table_filtered.qzv \
--verbose

# View on https://view.qiime2.org

#---------------------------------------------------------#

# Perform alignment on rep seqs

mkdir ~/gnoto/tree

qiime alignment mafft \
--i-sequences ~/gnoto/seqs/dada2out/gnoto_rep_seqs_filtered.qza \
--p-n-threads 5 \
--o-alignment ~/gnoto/tree/aligned_gnoto_rep_seqs.qza \
--verbose

#--------------------------------------------------#

# Mask highly variable regions of the alignment

qiime alignment mask \
--i-alignment ~/gnoto/tree/aligned_gnoto_rep_seqs.qza \
--o-masked-alignment ~/gnoto/tree/masked_aligned_gnoto_rep_seqs.qza \
--verbose

#--------------------------------------------------#

# Generate phylogenetic tree

qiime phylogeny fasttree \
--i-alignment ~/gnoto/tree/masked_aligned_gnoto_rep_seqs.qza \
--p-n-threads 1 \
--o-tree ~/gnoto/tree/gnoto_unrooted_tree.qza \
--verbose

#--------------------------------------------------#

# Apply mid-point rooting to tree

qiime phylogeny midpoint-root \
--i-tree ~/gnoto/tree/gnoto_unrooted_tree.qza \
--o-rooted-tree ~/gnoto/tree/gnoto_rooted_tree.qza \
--verbose

#--------------------------------------------------#

# Generate files for R

mkdir ~/gnoto/output

qiime tools export \
--input-path ~/gnoto/tree/gnoto_unrooted_tree.qza \
--output-path ~/gnoto/output/

qiime tools export \
--input-path ~/gnoto/taxonomy/gnoto_taxonomy.qza \
--output-path ~/gnoto/output/

qiime tools export \
--input-path ~/gnoto/seqs/dada2out/gnoto_table_filtered.qza \
--output-path ~/gnoto/output/

biom convert -i ~/gnoto/output/feature-table.biom \
--to-tsv \
-o ~/gnoto/output/gnoto_table.tsv

#--------------------------------------------------#

# Tidy up output files

# Remove header from feature table
sed -i '' '1d' ~/gnoto/output/gnoto_table.tsv

# Remove header from taxonomy table
sed -i '' '1,2d' ~/gnoto/output/gnoto_tax.tsv

# Remove confidence values from taxonomy table
cut -f 1,2 ~/gnoto/output/gnoto_tax.tsv > ~/gnoto/output/gnoto_tax.tsv

# Replace semi-colons in taxonomy table with tabs
# sed -i '' 's/;/<control><V><tab>/g' ~/gnoto/output/gnoto_tax.tsv

#-------------------- The end ---------------------#
