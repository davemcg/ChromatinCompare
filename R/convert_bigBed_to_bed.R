#' Wrapper to converted bigBed to bed, using UCSC bedToBigBed tool
#' 
#' @param A Path to the bigBed file
#' @details You must have the UCSC bedToBigBed tool installed. UCSC provides pre-built binaries here: http://hgdownload.soe.ucsc.edu/admin/exe/. \code{convert_bigBed_to_Bed} will return the path to the converted bed file. 
#' @examples bed_file_path <- convert_bigBed_to_bed('~/mm10_hela_chromhmm.bigBed')
#' 

convert_bigBed_to_bed <- function(bigBed){
  output_bed <- paste0(strsplit(bigBed,'\\.')[[1]][1], '.bed')
  
  if (!file.exists(output_bed)) {
    system(paste0('bigBedToBed ', bigBed, ' ', output_bed))
    system(paste0('bgzip ', output_bed))
  }
  
  return(paste0(output_bed,'.gz'))
}