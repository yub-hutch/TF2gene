#' Read ATAC-seq consensus peak
#'
#' This function reads a BED file and formats the consensus peaks. It only retains peaks in chromosome 1-22, X, & Y.
#'
#' @param fbed A character string specifying the path to the BED file.
#'
#' @return A tibble with four columns:
#' \describe{
#'   \item{chr}{Chromosome identifier.}
#'   \item{start}{Start position of the peak (0-based, included).}
#'   \item{end}{End position of the peak (0-based, excluded).}
#'   \item{openness}{Openness score of the peak.}
#' }
#' @export
#'
read_consensus_peaks <- function(fbed) {
  peak = dplyr::as_tibble(data.table::fread(fbed, header = F))[, c(1, 2, 3, 5)]
  names(peak) = c('chr', 'start', 'end', 'openness')
  return(peak)
}
