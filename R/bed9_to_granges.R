#' Converts imported bed file to granges object
#'
#' @param bed_dataframe Imported bed data.frame() object (from \code{bed9_importer})
#' @param genome Which genome is this? Defaults to mm10. Other options are FILL IN!!
#' @details \code{bed9_to_granges} Converts an imported bed file to a granges
#' object for annotation.
bed9_to_granges <- function(bed_dataframe, genome='mm10'){
  x <- bed_dataframe
  bed9_ranges <- GenomicRanges::GRanges(seqnames=x$chr, ranges = IRanges(start=as.numeric(x$start),end=as.numeric(x$end)), name=x$name, score = x$score, strand=(rep('*',length(x$chr))), thickStart = x$thickStart, thickEnd = x$thickEnd, itemRgb = x$itemRgb)
  GenomeInfoDb::genome(bed9_ranges) <- genome
  GenomeInfoDb::seqlevels(bed9_ranges) <- sort(GenomeInfoDb::seqlevels(bed9_ranges))
  return(bed9_ranges)
}
