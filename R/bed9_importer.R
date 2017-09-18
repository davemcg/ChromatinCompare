#' Imports bed9 (9 column) files into R
#'
#' @param bed Path to the bed file
#' @param uncompress_keep If the bed file is compressed, should the uncompressed file be
#' retained or removed? Defaults to TRUE.
#' @param skip How many lines to skip? Defaults to zero. You must look at the
#' first few lines of the bed file to determine if non bed entries are present.
#' @details \code{bed9_importer} imports a 9 column bed file and retains proper
#' encoding (character, integer) of the columns.
#' @example
#' bed9_importer('~/your_bed_file.bed.gz', skip=1 )
#' bed9_importer
bed9_importer <- function(bed, uncompress_keep = TRUE, skip=0){
  not_gzip <- strsplit(bed,'\\.gz')[[1]][1]
  if (file.exists(not_gzip)){
    output <- readr::read_tsv(not_gzip, col_names = F, col_types = list('c','i','i','c','i','c','i','i','c'), skip = skip)
  } else {
    output <- readr::read_tsv(bed, col_names = F, col_types = list('c','i','i','c','i','c','i','i','c'), skip = skip)
  }
  if (!uncompress_keep) {
    system(paste0('rm ', not_gzip))
  }
  colnames(output) <- c('chr','start','end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb')
  return(output)
}
