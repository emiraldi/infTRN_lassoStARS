# build_prior_utils.R
# Functions used for building and processing prior networks.


# peaks_region_2_name
# Convert dataframe defining genomic regions into character
# vector of region names
# Inputs:
# df_peak = 3 column dataframe of the form [Chr Start End]
# delimiter = delimiter for elements in character vector (default: '_')
# Output:
# name_peak = character vector of peak names
peaks_region_2_name <- function(
    df_peak,
    delimiter='_'
    ){
    name_peak = paste(df_peak[,1], df_peak[,2], df_peak[,3], sep=delimiter)
    return(name_peak)
}


# net_sparse_2_full
# Convert sparse network to full matrix format
# Network should be of the form: [Node1, Node2, Weight, ...]
# Only the first three columns will be used.
# Inputs:
# net_sp = sparse network file name or sparse file
# Outputs:
# net_full = full matrix network, rows=Node2, columns=Node1
net_sparse_2_full <- function(
    net_sp
    ){
    
    library('Matrix')
    
    # load sparse network if necessary
    if (is.character(net_sp)){
        print(paste('Load sparse network:',net_sp))
        net_sp <- read.delim(net_sp, header=TRUE, sep='\t')
    }
    net_sp <- net_sp[,1:3]
    colnames(net_sp) <- c('Node1', 'Node2', 'Weight')
    
    # unique TFs & Targets
    uniq_node1 <- unique(net_sp$Node1)
    uniq_node2 <- unique(net_sp$Node2)
    n_uniq_node1 <- length(uniq_node1)
    n_uniq_node2 <- length(uniq_node2)
    
    # assign unique nodes to indices
    df_node1 <- data.frame(idx=1:n_uniq_node1, row.names=uniq_node1)
    df_node2 <- data.frame(idx=1:n_uniq_node2, row.names=uniq_node2)
    
    # create sparse network using indices
    idx_node1 <- df_node1[as.character(net_sp$Node1),1]
    idx_node2 <- df_node2[as.character(net_sp$Node2),1]
    net_sp_idx <- data.frame(Node1=idx_node1, Node2=idx_node2)
    net_sp_idx <- cbind(net_sp_idx, net_sp[,3])
    colnames(net_sp_idx) <- c('Node1', 'Node2', 'Weight')

    # full network
    net_full <- sparseMatrix(i=net_sp_idx$Node2, j=net_sp_idx$Node1, x=net_sp_idx$Weight)
    net_full <- as.matrix(net_full)
    rownames(net_full) <- as.character(uniq_node2)
    colnames(net_full) <- as.character(uniq_node1)
    
    return(net_full)
}


# net_full_2_sparse
# Convert network in full matrix format to sparse format
# Sparse network will have the form [TF, Target, Weight]
# where TF and Target are the columns and rows of the full
# matrix, respectively, and Weight are the non-zero matrix entries.
# Input:
# net_full <- full network matrix, tab delimited, rows=Targets, cols=TFs
# Output:
# net_sp <- sparse network, 3 column format: [TF, Target, Weight]
net_full_2_sparse <- function(
    net_full
    ){
    
    # sparse matrix
    net_sp <- as.data.frame(as.table(as.matrix(net_full)))
    net_sp <- data.frame(TF=net_sp[,2], Target=net_sp[,1], Weight=net_sp[,3])
    
    # delete zero entries
    idx_zero <- which(net_sp$Weight==0)
    if (length(idx_zero) > 0){
        net_sp <- net_sp[-idx_zero,]
    }
    
    return(net_sp)
}


# net_quant_2_binary
# Convert sparse network to binary network
# All interaction weights will be set to 1
# Network should be of the form: [TF, Target, Weight, ...]
# Only the first three columns will be used.
# Inputs:
# net_sp = sparse network, 3 column format: [TF, Target, Weight]
# Outputs:
# net_sp_b = sparse network, all weights equal to 1
net_quant_2_binary <- function(
    net_sp
    ){

    # load sparse network if necessary
    if (is.character(net_sp)){
        print(paste('Load sparse network:',net_sp))
        net_sp <- read.delim(net_sp, sep='\t', header=TRUE)
    }
    net_sp <- net_sp[,1:3]
    
    # Set weights to 1
    net_sp_b <- data.frame(TF=net_sp[,1], Target=net_sp[,2], Weight=1)
    
    return(net_sp_b)
}


# net_sum
# Sum of two full networks
# Dimensions do not have to be the same
# Inputs:
# net1 = network 1, full matrix format
# net2 = network 2, full matrix format
# Outputs:
# net_out = sum of input networks
net_sum <- function(
    net1,
    net2
    ){
    
    # unique rownames and colnames
    row_names <- sort(unique(c(rownames(net1), rownames(net2))))
    n_row <- length(row_names)
    col_names <- sort(unique(c(colnames(net1), colnames(net2))))
    n_col <- length(col_names)

    # initialize sum matrix
    net1_sum <- matrix(0,n_row,n_col)
    rownames(net1_sum) <- row_names
    colnames(net1_sum) <- col_names
    net2_sum <- matrix(0,n_row,n_col)
    rownames(net2_sum) <- row_names
    colnames(net2_sum) <- col_names

    # sum inputs
    net1_sum[rownames(net1),colnames(net1)] <- net1
    net2_sum[rownames(net2),colnames(net2)] <- net2
    net_out <- net1_sum + net2_sum
    
    return(net_out)
}


# prior_proximal
# Build prior based on motif occurrences proximal to a region
# surrounding genomic features.
# Inputs:
# bed_motif = bed file/dataframe of motif locations, 4-column format: [Chr, Start, End, Motif]
# bed_feature = bed file/dataframe of genomic features to map TFs to: 4-column format: [Chr, Start, End, Gene]
# tf_motif = file/dataframe mapping motifs to TF names, 2-column format: [Motif, TF]
# window_feature = window (# kb) around bed features to map to TFs (default: 10000)
# bed_active = bed file/dataframe of active genomic regions (eg, active histone markers) (default: NULL)
# Output:
# prior = 3-column prior network [TF Target Weight]
prior_proximal <- function(
    bed_motif,
    bed_feature,
    tf_motif,
    window_feature=10000,
    bed_active=NULL
    ){
    
    # load motif bed file
    if (is.character(bed_motif)){
        print(paste('Load motif bed file:',bed_motif))
        bed_motif <- read.delim(bed_motif, header=FALSE, sep='\t')
    }
    bed_motif <- bed_motif[,1:4]
    colnames(bed_motif) <- c('Chr', 'Start', 'End', 'Motif')
    
    # load feature bed file
    if (is.character(bed_feature)){
        print(paste('Load feature bed file:',bed_feature))
        bed_feature <- read.delim(bed_feature, header=FALSE, sep='\t')
    }
    bed_feature <- bed_feature[,1:4]
    
    # load TF to motif map if necessary
    if (is.character(tf_motif)){
        print(paste('Load TF to motif info:',tf_motif))
        tf_motif <- read.delim(tf_motif, header=FALSE, sep='\t')
    }
    tf_motif <- tf_motif[,1:2]
    colnames(tf_motif) <- c('Motif','TF')
    tf_motif <- data.frame(TF=tf_motif$TF, Motif=tf_motif$Motif)
    
    # prune bed features if necessary
    if (!is.null(bed_active)){
        # load active bed
        if (is.character(bed_active)){
            print(paste('Load active bed file:',bed_active))
            bed_active <- read.delim(bed_active, header=FALSE, sep='\t')
        }
        print(paste('Prune feature bed using active bed'))
        bed_feature <- prune_bed_feature(bed_feature=bed_feature, bed_keep_these=bed_active)
    }
    colnames(bed_feature) <- c('Chr', 'Start', 'End', 'Gene')
    
    # create temporary working directory
    dir_wd <- basename(tempfile(pattern='tmp_dir_'))
    dir.create(dir_wd, showWarnings=TRUE)
    
    # create motif bed file
    print('Create motif bed file')
    write.table(merge(bed_motif, tf_motif, by='Motif')[,2:5], file.path(dir_wd,'tmp_bed_motif.bed'), 
                    quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')
    
    # sort motif bed file
    system(paste0('bedtools sort -i ', file.path(dir_wd,'tmp_bed_motif.bed'),' > ', 
                    file.path(dir_wd,'tmp_bed_motif_sorted.bed')))
    
    # save bed features to working
    write.table(bed_feature, file.path(dir_wd,'tmp_bed_feature.bed'), quote=FALSE, sep='\t',
                    row.names=FALSE, col.names=FALSE)
    
    # sort bed features file
    system(paste0('bedtools sort -i ', file.path(dir_wd,'tmp_bed_feature.bed'),' > ', 
                    file.path(dir_wd,'tmp_bed_feature_sorted.bed')))
    
    # intersection of gene features and motifs
    print('Window motifs with gene features')
    cmd_window <- paste0("bedtools window -w ",window_feature," -a ",
                            dir_wd,"/tmp_bed_feature_sorted.bed -b ",
                            dir_wd,"/tmp_bed_motif_sorted.bed | awk -v OFS='\t' '{print $8,$4}' > ",
                            dir_wd,"/tmp_window_interactions.bed")
    system(cmd_window)
    
    # load interactions
    df_int <- read.delim(file.path(dir_wd,'tmp_window_interactions.bed'), header=FALSE, sep='\t')
    colnames(df_int) <- c('TF','Target')
    
    # delete temporary working directory
    system(paste0('rm -r ',dir_wd))
    
    # count number of motifs in vicinity of each gene
    print('Count TF-gene interactions')
    name_int <- paste(df_int[,1], df_int[,2], sep='_')
    idx_unique <- which(!duplicated(name_int))
    unique_int <- df_int[idx_unique,]
    name_unique_int <- paste(unique_int[,1], unique_int[,2], sep='_')
    df_name_int <- as.data.frame(table(name_int))
    df_name_int <- data.frame(Freq=df_name_int$Freq, row.names=df_name_int$name_int)
    
    # build prior
    print('Build prior')
    prior <- data.frame(TF=unique_int$TF, Target=unique_int$Target, Weight=df_name_int[name_unique_int,1])
    name_prior <- paste(prior$TF, prior$Target, sep='_')
    sort_prior <- sort(name_prior, index.return=TRUE)
    prior <- prior[sort_prior$ix,]
    rownames(prior) <- 1:dim(prior)[1]
    
    return(prior)
}


# prune_bed_feature
# Returns features specified by a bed file that overlap another bed file.
# Inputs:
# bed_feature = bed file of genomic features to prune
# bed_keep_these = keep the feature in bed_feature that overlap with these
# keep_unique = keep unique pruned bed file features (default: TRUE)
# Outputs
# bed_pruned = bed file of genomic features overlapped desired region
prune_bed_feature <- function(
    bed_feature,
    bed_keep_these,
    keep_unique=TRUE
    ){

    # create temporary working directory
    dir_wd <- basename(tempfile(pattern='tmp_dir_'))
    dir.create(dir_wd, showWarnings=TRUE)
    
    # load bed files
    if (is.character(bed_feature)){
        bed_feature <- read.delim(bed_feature, header=FALSE, sep='\t')
    }
    if (is.character(bed_keep_these)){
        bed_keep_these <- read.delim(bed_keep_these, header=FALSE, sep='\t')
    }
    
    # number of input bed columns
    n_col_bed_feature <- dim(bed_feature)[2]
    n_col_bed_keep_these <- dim(bed_keep_these)[2]
    
    # write bed files to temporary working directory
    write.table(bed_feature, file.path(dir_wd, 'bed_feature.bed'), quote=FALSE,
                    sep='\t', col.names=FALSE, row.names=FALSE)
    write.table(bed_keep_these, file.path(dir_wd, 'bed_keep_these.bed'), quote=FALSE,
                    sep='\t', col.names=FALSE, row.names=FALSE)
    
    # sort bed files
    system(paste0('bedtools sort -i ',file.path(dir_wd, 'bed_feature.bed'),' > ', 
                        file.path(dir_wd, 'bed_feature_sorted.bed')))
    system(paste0('bedtools sort -i ',file.path(dir_wd, 'bed_keep_these.bed'),' > ', 
                        file.path(dir_wd, 'bed_keep_these_sorted.bed')))
    
    # prune bed features
    cmd_window <- paste0('bedtools window -w 0 -a ', file.path(dir_wd,'bed_keep_these_sorted.bed'), ' -b ',
                            file.path(dir_wd, 'bed_feature.bed'),' > ', file.path(dir_wd,'bed_window.bed'))
    system(cmd_window)
    n_line_window <- length(readLines(file.path(dir_wd,'bed_window.bed')))
    if (n_line_window > 0){
        bed_pruned <- read.delim(file.path(dir_wd,'bed_window.bed'), header=FALSE, sep='\t')
        bed_pruned <- bed_pruned[,(n_col_bed_keep_these+1):(n_col_bed_keep_these+n_col_bed_feature)]
        # keep unique features
        if (keep_unique){
            name_feat <- NULL
            for (ix in 1:n_col_bed_feature){
                name_feat <- paste(name_feat, bed_pruned[,ix], sep='_')
            }
            idx_keep <- which(!duplicated(name_feat))
            bed_pruned <- bed_pruned[idx_keep,]
        }
    } else {
        bed_pruned <- NULL
    }
    
	# remove temporary working directory
    system(paste0('rm -r ',dir_wd))
    
    return(bed_pruned)
}


# save_data_matrix
# Save gene expression matrix.
# txt file format, tab-delimited, rows=genes, cols=cells
# Inputs:
# counts_mat = gene expression matrix
# filename = output file name
# out_digits = maximum number of significant digits to keep (default: 5)
save_data_matrix <- function(
    counts_mat,
    filename,
    out_digits=5
    ){
    write.table(signif(counts_mat, digits=out_digits), filename, quote=FALSE, row.names=TRUE, col.names=NA, sep='\t')
}


# load_data_matrix
# Load gene expression matrix.
# txt file format, tab-delimited, rows=genes, cols=cells
# header are unique cell names, first column is gene names
# Inputs:
# filename = filename string
# replace_na = replace NAs with zeros (default: FALSE)
# Output:
# counts_mat = gene expression matrix
load_data_matrix <- function(
    filename,
    replace_na=FALSE
    ){

    counts_mat <- read.delim(filename, header=TRUE, row.names=1, check.names=FALSE)
    if (replace_na){
        counts_mat[is.na(counts_mat)] <- 0
    }
    counts_mat <- as.matrix(counts_mat)

    return(counts_mat)
}


# net_prior_merge
# Create merged prior.
# Inputs:
# net_full = full network, format: [Targets X TFs]
# Output:
# net_out = list containing full network and degenerate TF info
net_prior_merge <- function(
    net_full
    ){

    # load network
    if (is.character(net_full)){
        print(paste('Load network:',net_full))
        net_full <- load_data_matrix(net_full)
    }
    
	# unique string for each TF
	n_tf <- dim(net_full)[2]
	str_tf <- NULL
	for (ix in 1:n_tf){
		curr_str <- paste(net_full[,ix], collapse='|')
		str_tf <- c(str_tf, curr_str)
	}

	# group TFs with same target profile
	uniq_str_tf <- unique(str_tf)
	n_uniq_str <- length(uniq_str_tf)
	net_merged <- NULL
	all_name <- NULL
	merged_tfs <- NULL
	for (ix in 1:n_uniq_str){
		curr_str <- uniq_str_tf[ix]
		idx_tf <- which(str_tf==curr_str)
		curr_name <- paste(colnames(net_full)[idx_tf], collapse='_')
		all_name <- c(all_name,curr_name)
		net_merged <- cbind(net_merged, net_full[,idx_tf[1]])
		if (length(idx_tf)>1){
			curr_merged <- data.frame(name_merged=curr_name, 
			                            tfs_merged=paste(colnames(net_full)[idx_tf], collapse=', '))
			merged_tfs <- rbind(merged_tfs, curr_merged)
		}
	}
	colnames(net_merged) <- all_name
	rownames(net_merged) <- rownames(net_full)
	
	net_out <- list()
	net_out[['Network']] <- net_merged
	net_out[['MergedTFs']] <- merged_tfs
	return(net_out)
	
}
