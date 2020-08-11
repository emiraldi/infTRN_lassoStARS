# construct_atac_prior.R
# Generate prior gene regulatory network from: 
# 1. bed file of accessible chromatin regions.
# 2. TF motifs
# 3. bed file of genomic features (eg, TSS, gene body)
print('--------------------------------------')
rm(list=ls())
options(stringsAsFactors=FALSE)

library(motifmatchr)
library(GenomicRanges)
library(TFBSTools)
source('./utils_prior.R')
source('./utils_bedtools.R')

#================== INPUTS ===================

# output directory
dir_out <- './example_prior_outputs'

# prior name (used for output filenames)
name_prior <- 'prior_example_atac'

# bed file of ATAC-seq peaks
file_peaks <- './example_prior_inputs/example_atac_peaks.bed'

# genome build
genome <- 'hg19'

# motif directory
dir_motif <- './example_prior_inputs/example_motifs'

# file mapping motifs to genes
# 2-column: [motif name, gene name]
file_motif_info <- './example_prior_inputs/tbl_motif_2_gene.tsv'

# motif scanning p-value cutoff
pval_cutoff <- 1E-5

# genomic feature file (e.g., gene body, TSS)
# 4-column: [chr, start, end, feature name]
file_features <- './example_prior_inputs/example_gene_body.bed'

# window within genomic features to search (eg, +/-10kb gene body)
window_feature <- 10000

#=============================================

# create output directory
print(paste('Create output directory:', dir_out))
dir.create(dir_out, recursive=TRUE, showWarnings=FALSE)

# load peaks
print(paste('Load peaks:', file_peaks))
peaks <- read.delim(file_peaks, header=FALSE, sep='\t')
peaks <- peaks[,1:3]

# load motifs
print(paste('Load motifs in directory:', dir_motif))
motifs <- list()
for (ix in list.files(dir_motif,full.names=TRUE)){
    curr_motif <- as.matrix(read.delim(ix, header=FALSE, sep='\t'))
    curr_name <- tools::file_path_sans_ext(basename(ix))
    rownames(curr_motif) <- c('A','C','G','T')
    motifs[[curr_name]] <- PFMatrix(ID=curr_name, name=curr_name, bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                                        profileMatrix=100*curr_motif)
}
motifs <- do.call(PFMatrixList, motifs)

# load motif to gene mapping
print(paste('Load motif info:', file_motif_info))
motif_info <- read.delim(file_motif_info, header=FALSE, sep='\t')
colnames(motif_info) <- c('Motif','TF')

# load genomic features
print(paste('Load genomic features:', file_features))
features <- read.delim(file_features, header=FALSE, sep='\t')

# filter peaks for those within window of genomic features
print('Filter peaks')
features_expand <- features[,1:3]
features_expand[,2] <- pmax(features_expand[,2]-window_feature,0)
features_expand[,3] <- features_expand[,3] + window_feature
peaks_filtered <- bedtools_intersect(peaks, features_expand)
peaks_filtered <- bedtools_merge(peaks_filtered)
peaks_range <- GRanges(seqnames=peaks_filtered[,1], 
                        ranges=IRanges(start=peaks_filtered[,2], end=peaks_filtered[,3]))

# motif scanning
print('Scan peaks for motifs')
motif_scan <- matchMotifs(motifs, peaks_range, genome=genome, p.cutoff=pval_cutoff, out='positions')
motif_scan <- as.data.frame(motif_scan)[,c('seqnames', 'start', 'end', 'group_name')]
colnames(motif_scan) <- c('Chr','Start','End','Motif')

# construct prior - quantitative, sparse format
print('Construct prior')
prior_q_sp <- prior_proximal(bed_motif=motif_scan, bed_feature=features,
                                tf_motif=motif_info, window_feature=window_feature)

print('Save priors')
# quantitative - sparse
file_out <- file.path(dir_out, paste0(name_prior,'_q_sp.tsv'))
write.table(prior_q_sp, file_out, quote=FALSE, sep='\t', col.names=TRUE, row.names=FALSE)
# quantitative - full
prior_q_full <- net_sparse_2_full(prior_q_sp)
file_out <- file.path(dir_out, paste0(name_prior,'_q.tsv'))
save_data_matrix(prior_q_full, file_out)
# binary - sparse
prior_b_sp <- net_quant_2_binary(prior_q_sp)
file_out <- file.path(dir_out, paste0(name_prior,'_b_sp.tsv'))
write.table(prior_b_sp, file_out, quote=FALSE, sep='\t', col.names=TRUE, row.names=FALSE)
# binary - full
prior_b_full <- net_sparse_2_full(prior_b_sp)
file_out <- file.path(dir_out, paste0(name_prior,'_b.tsv'))
save_data_matrix(prior_b_full, file_out)

print('Save merged priors')
# Given degeneracy of TF motifs, many priors based on ATAC-seq data 
# contain TFs with identical target genes and interaction strengths.
# We merge these TFs to create meta-TFs.
# merged quantitative - full
prior_merge_q_full <- net_prior_merge(prior_q_full)
if ('MergedTFs' %in% names(prior_merge_q_full)){
    file_out <- file.path(dir_out, paste0(name_prior,'_q_merged.tsv'))
    save_data_matrix(prior_merge_q_full[['Network']], file_out)
    # merged quantitative - sparse
    prior_merge_q_sp <- net_full_2_sparse(prior_merge_q_full[['Network']])
    file_out <- file.path(dir_out, paste0(name_prior,'_q_merged_sp.tsv'))
    write.table(prior_merge_q_sp, file_out, quote=FALSE, sep='\t', col.names=TRUE, row.names=FALSE)
    # merged quantitative - merged TF info
    file_out <- file.path(dir_out, paste0(name_prior,'_q_mergedTfs.txt'))
    write.table(prior_merge_q_full[['MergedTFs']], file_out, quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')
}
# merged binary - sparse
prior_merge_b_full <- net_prior_merge(prior_b_full)
if ('MergedTFs' %in% names(prior_merge_b_full)){
    file_out <- file.path(dir_out, paste0(name_prior,'_b_merged.tsv'))
    save_data_matrix(prior_merge_b_full[['Network']], file_out)
    # merged binary - full
    prior_merge_b_sp <- net_full_2_sparse(prior_merge_b_full[['Network']])
    file_out <- file.path(dir_out, paste0(name_prior,'_b_merged_sp.tsv'))
    write.table(prior_merge_b_sp, file_out, quote=FALSE, sep='\t', col.names=TRUE, row.names=FALSE)
    # merged quantitative - merged TF info
    file_out <- file.path(dir_out, paste0(name_prior,'_b_mergedTfs.txt'))
    write.table(prior_merge_b_full[['MergedTFs']], file_out, quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')
}

#=============================================
print('--------------------------------------')
print('Done!')