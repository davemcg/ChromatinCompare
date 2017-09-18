#' Load species and build specific annotations
#'
#' @param bed Species and genome build
#' @param type Gene body is the default. You can also load exon only, promoter
#' only, intron only.

#' @details \code{bed9_importer} imports a 9 column bed file and retains proper
#' encoding (character, integer) of the columns.
#' @example
#' bed9_importer('~/your_bed_file.bed.gz', skip=1 )
#' bed9_importer
annots = c('mm10_basicgenes')
annotations = build_annotations(genome = 'mm10', annotations = annots)
seqlevelsStyle(annotations) <- "UCSC"
seqlevels(annotations) <- sort(seqlevels(annotations))
# remove annotations without a gene symbol and LOC* and LINC*
annotations <- annotations %>% data.frame() %>% filter(!is.na(symbol)) %>% filter(!grepl('^LOC', symbol)) %>% filter(!grepl('^LINC',symbol))
# intergenic
annots_intergenic <- c('mm10_genes_intergenic')
intergenic_annotations = build_annotations(genome = 'mm10', annotations = annots_intergenic)
seqlevelsStyle(intergenic_annotations) <- "UCSC"
seqlevels(intergenic_annotations) <- sort(seqlevels(intergenic_annotations))

# merge annotations with intergenic
annotations <- rbind(annotations, intergenic_annotations %>% data.frame())
# rebuild granges
annotations <- GRanges(seqnames=annotations$seqnames, ranges=IRanges(start=annotations$start, end = annotations$end), strand = annotations$strand , id = annotations$id, gene_id = annotations$id, symbol = annotations$symbol, type =annotations$type)
seqlevels(annotations) <- sort(seqlevels(annotations))
