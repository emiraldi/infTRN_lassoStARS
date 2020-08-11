# utils_bedtools.R
# Functions to perform bedtools operations in R.


# bedtools_sort
# Sort bed file.
# Inputs:
# bed_input = bed file to sort: dataframe with first 3 columns [chr start end]
# Output:
# bed_sorted = sorted bed file: dataframe with first 3 columns [chr start end]
bedtools_sort <- function(
    bed_input
    ){
    
    # original column names
    name_col <- colnames(bed_input)
    
    # create temporary working directory
    dir_tmp <- tempfile(pattern='dir', fileext='')
    dir.create(dir_tmp)
    
    # make sure writing integers to file
    bed_input[,2] <- as.integer(bed_input[,2])
    bed_input[,3] <- as.integer(bed_input[,3])
    
    # write bed to file
    file_unsorted <- file.path(dir_tmp,'unsorted.bed')
    write.table(bed_input, file_unsorted, quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
    
    # sort bed
    file_sorted <- file.path(dir_tmp,'sorted.bed')
    cmd_sort <- paste0('bedtools sort -i ',file_unsorted,' > ', file_sorted)
    system(cmd_sort)
    
    # load sorted bed
    bed_sorted <- read.delim(file_sorted, header=FALSE, sep='\t')
    colnames(bed_sorted) <- name_col
    
    # remove temporary working directory
    system(paste0('rm -r ',dir_tmp))
    
    return(bed_sorted)
}


# bedtools_merge
# Merge overlapping regions of a bed file.
# Inputs:
# bed_input = bed file to merge: dataframe with 3 columns [chr start end]
# is_sorted = is the input bed file sorted? (default: FALSE)
# Output:
# bed_merged = merged bed file: dataframe with 3 columns [chr start end]
bedtools_merge <- function(
    bed_input,
    is_sorted=FALSE
    ){
    
    # return error if input is not 3 columns
    n_col <- dim(bed_input)[2]
    if (n_col != 3){
        stop('Input bed dataframe must have 3 columns')
    }
    
    # original column names
    name_col <- colnames(bed_input)
    
    # create temporary working directory
    dir_tmp <- tempfile(pattern='dir', fileext='')
    dir.create(dir_tmp)
    
    # make sure writing integers to file
    bed_input[,2] <- as.integer(bed_input[,2])
    bed_input[,3] <- as.integer(bed_input[,3])
    
    # write input bed to file
    file_sorted <- file.path(dir_tmp,'sorted.bed')
    if (is_sorted){
        write.table(bed_input, file_sorted, quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
    } else {
        # sort input bed file
        file_unsorted <- file.path(dir_tmp,'unsorted.bed')
        write.table(bed_input, file_unsorted, quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
        cmd_sort <- paste0('bedtools sort -i ',file_unsorted,' > ', file_sorted)
        system(cmd_sort)
    }
    
    # merge bed
    file_merged <- file.path(dir_tmp,'merged.bed')
    cmd_merge <- paste0('bedtools merge -i ',file_sorted,' > ', file_merged)
    system(cmd_merge)
    
    # load merged bed
    bed_merged <- read.delim(file_merged, header=FALSE, sep='\t')
    colnames(bed_merged) <- name_col
    
    # remove temporary working directory
    system(paste0('rm -r ',dir_tmp))
    
    return(bed_merged)
}


# bedtools_window
# Search for overlapping features between two bed files
# Inputs:
# bed_A = bed file A: dataframe with first 3 columns [chr start end]
# bed_B = bed file B: dataframe with first 3 columns [chr start end]
# window_A = window to add upstream and downstream of A (default: 0)
# is_sorted = are the input bed files sorted? (default: FALSE)
# keep_cols = maintain all columns of input bed instead of 1:3 (default: FALSE)
# Output:
# bed_window = features in A that overlap those in B, 6 column dataframe
bedtools_window <- function(
    bed_A,
    bed_B,
    window_A=0,
    is_sorted=FALSE,
    keep_cols=FALSE
    ){
    
    # create temporary working directory
    dir_tmp <- tempfile(pattern='dir', fileext='')
    dir.create(dir_tmp)
    
    # Use only first 3 columns
    if (!keep_cols){
        bed_A <- bed_A[,1:3]
        bed_B <- bed_B[,1:3]
    }
    
    # make sure writing integers to file
    bed_A[,2] <- as.integer(bed_A[,2])
    bed_A[,3] <- as.integer(bed_A[,3])
    bed_B[,2] <- as.integer(bed_B[,2])
    bed_B[,3] <- as.integer(bed_B[,3])
    
    # sort if necessary
    if (!is_sorted){
        bed_A <- bedtools_sort(bed_A)
        bed_B <- bedtools_sort(bed_B)
        bed_A[,2] <- as.integer(bed_A[,2])
        bed_A[,3] <- as.integer(bed_A[,3])
        bed_B[,2] <- as.integer(bed_B[,2])
        bed_B[,3] <- as.integer(bed_B[,3])
    }
    
    # write input bed to file
    file_sorted_A <- file.path(dir_tmp,'sorted_A.bed')
    file_sorted_B <- file.path(dir_tmp,'sorted_B.bed')
    write.table(bed_A, file_sorted_A, quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
    write.table(bed_B, file_sorted_B, quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
    
    # window bed
    file_window <- file.path(dir_tmp,'window.bed')
    cmd_window <- paste0('bedtools window -w ',window_A,' -a ',file_sorted_A,' -b ',file_sorted_B,' > ', file_window)
    system(cmd_window)
    
    # load window bed
    bed_window <- read.delim(file_window, header=FALSE, sep='\t')
    
    # remove temporary working directory
    system(paste0('rm -r ',dir_tmp))
    
    return(bed_window)
}


# bedtools_intersect
# Search for intersecting features between two bed files
# Inputs:
# bed_A = bed file A: dataframe with first 3 columns [chr start end]
# bed_B = bed file B: dataframe with first 3 columns [chr start end]
# is_sorted = are the input bed files sorted? (default: FALSE)
# wa = Write the original entry in A for each overlap (default: FALSE)
# wb = Write the original entry in B for each overlap (default: FALSE)
# Output:
# bed_intersect = intersection of bed A and bed B, 3 or 6 column dataframe
bedtools_intersect <- function(
    bed_A,
    bed_B,
    is_sorted=FALSE,
    wa=FALSE,
    wb=FALSE
    ){
    
    # create temporary working directory
    dir_tmp <- tempfile(pattern='dir', fileext='')
    dir.create(dir_tmp)
    
    # Use only first 3 columns
    bed_A <- bed_A[,1:3]
    bed_B <- bed_B[,1:3]
    
    # make sure writing integers to file
    bed_A[,2] <- as.integer(bed_A[,2])
    bed_A[,3] <- as.integer(bed_A[,3])
    bed_B[,2] <- as.integer(bed_B[,2])
    bed_B[,3] <- as.integer(bed_B[,3])
    
    # sort if necessary
    if (!is_sorted){
        bed_A <- bedtools_sort(bed_A)
        bed_B <- bedtools_sort(bed_B)
        bed_A[,2] <- as.integer(bed_A[,2])
        bed_A[,3] <- as.integer(bed_A[,3])
        bed_B[,2] <- as.integer(bed_B[,2])
        bed_B[,3] <- as.integer(bed_B[,3])
    }
    
    # write input bed to file
    file_sorted_A <- file.path(dir_tmp,'sorted_A.bed')
    file_sorted_B <- file.path(dir_tmp,'sorted_B.bed')
    write.table(bed_A, file_sorted_A, quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
    write.table(bed_B, file_sorted_B, quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
    
    # intersect bed
    file_intersect <- file.path(dir_tmp,'intersect.bed')
    cmd_intersect <- paste0('bedtools intersect -a ',file_sorted_A,' -b ',file_sorted_B)
    if (wa){
        cmd_intersect <- paste0(cmd_intersect,' -wa')
    }
    if (wb){
        cmd_intersect <- paste0(cmd_intersect,' -wb')
    }
    cmd_intersect <- paste0(cmd_intersect,' > ', file_intersect)
    system(cmd_intersect)
    
    # load intersect bed
    bed_intersect <- read.delim(file_intersect, header=FALSE, sep='\t')
    
    # remove temporary working directory
    system(paste0('rm -r ',dir_tmp))
    
    return(bed_intersect)
}


# bedtools_set_size
# Set all regions of bed file to specified size.
# Inputs:
# bed_input = bed file: dataframe with first 3 columns [chr start end]
# fix_size = fixed size of output regions (default: 500)
# fix_point = coordinate to fix (default: 'center')
#   options: 'center', 'start', 'end'
# Output:
# bed_fixed = bed file: dataframe with first 3 columns [chr start end]
bedtools_set_size <- function(
    bed_input,
    fix_size=500,
    fix_point='center'
    ){
    
    library('rtracklayer')
    library('GenomicRanges')
    
    # original column names
    name_col <- colnames(bed_input)
    
    # create temporary working directory
    dir_tmp <- tempfile(pattern='dir', fileext='')
    dir.create(dir_tmp)
    
    # make sure writing integers to file
    bed_input[,2] <- as.integer(bed_input[,2])
    bed_input[,3] <- as.integer(bed_input[,3])
    
    # write bed to file
    file_inputbed <- file.path(dir_tmp,'inputbed.bed')
    write.table(bed_input, file_inputbed, quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
    
    # fix size of regions
    user_ranges <- import.bed(file_inputbed)
    resize_ranges <- resize(user_ranges, width=fix_size, fix=fix_point)
    file_outputbed <- file.path(dir_tmp,'outputbed.bed')
    export.bed(resize_ranges, file_outputbed, format='bed')
    
    # load fixed bed
    bed_fixed <- read.delim(file_outputbed, header=FALSE, sep='\t')
    bed_fixed <- bed_fixed[,1:3]
    colnames(bed_fixed) <- name_col
    
    # remove temporary working directory
    system(paste0('rm -r ',dir_tmp))
    
    return(bed_fixed)
}

