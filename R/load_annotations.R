#' Load species and build specific annotations
#'
#' @param bed Species and genome build
#' @param type Gene body is the default. You can also load exon only, promoter
#' only, intron only.
#' @details \code{load_annotations} Loads in species and build specific annotations
#' encoding (character, integer) of the columns.
#' @example
#' annotations <- load_annotations(build= 'mm10', type='gene_byd' )
#' annotations
load_annotations <- function(build='mm10', type='gene_body'){
  annots = c('mm10_basicgenes')
  annotations = annotatr::build_annotations(genome = 'mm10', annotations = annots)
  seqlevelsStyle(annotations) <- "UCSC"
  GenomeInfodb::seqlevels(annotations) <- sort(GenomeInfodb::seqlevels(annotations))
  # remove annotations without a gene symbol and LOC* and LINC*
  annotations <- annotations %>% data.frame() %>% filter(!is.na(symbol)) %>% filter(!grepl('^LOC', symbol)) %>% filter(!grepl('^LINC',symbol))
  # intergenic
  annots_intergenic <- c('mm10_genes_intergenic')
  intergenic_annotations = annotatr::build_annotations(genome = 'mm10', annotations = annots_intergenic)
  GenomeInfodb::seqlevelsStyle(intergenic_annotations) <- "UCSC"
  GenomeInfodb::seqlevels(intergenic_annotations) <- sort(GenomeInfodb::seqlevels(intergenic_annotations))

  # merge annotations with intergenic
  annotations <- rbind(annotations, intergenic_annotations %>% data.frame())
  # rebuild granges
  annotations <- GenomicRanges::GRanges(seqnames=annotations$seqnames, ranges=IRanges(start=annotations$start, end = annotations$end), strand = annotations$strand , id = annotations$id, gene_id = annotations$id, symbol = annotations$symbol, type =annotations$type)
  GenomeInfodb::seqlevels(annotations) <- sort(GenomeInfodb::seqlevels(annotations))
  return(annotations)
}
