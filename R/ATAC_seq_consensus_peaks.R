#' Format ATAC-seq consensus peak
#'
#' This function reads a BED file and formats the consensus peaks. It only retains peaks in chromosome 1-22, X, & Y.
#'
#' @param fbed A character string specifying the path to the BED file.
#'
#' @return A tibble with four columns:
#' \describe{
#'   \item{chr}{Chromosome identifier.}
#'   \item{start}{Start position of the peak.}
#'   \item{end}{End position of the peak.}
#'   \item{openness}{Openness score of the peak.}
#' }
#' @export
#'
format_consensus_peaks <- function(fbed) {
  peak = dplyr::as_tibble(data.table::fread(fbed, header = F))[, c(1, 2, 3, 5)]
  names(peak) = c('chr', 'start', 'end', 'openness')

  peak = peak %>%
    dplyr::filter(chr %in% paste0('chr', c(as.character(1:22), 'X', 'Y'))) %>%
    dplyr::mutate(chr = sapply(strsplit(chr, 'chr'), tail, n = 1))

  return(peak)
}
