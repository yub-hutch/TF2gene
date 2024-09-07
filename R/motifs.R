#' Get Cluster-Buster Motif Information
#'
#' This function retrieves information about motifs from a specified directory.
#'
#' @param dir_motif A character string specifying the directory containing motif files.
#' @param ncores An integer specifying the number of cores to use for parallel processing.
#' @param verbose A logical value indicating whether to display progress. Default is TRUE.
#'
#' @return A tibble with the following columns:
#' \describe{
#'   \item{motif}{A character vector containing motif IDs.}
#'   \item{num_inner_motif}{An integer vector containing the number of inner motifs for each motif ID.}
#' }
#'
#' @export
get_cb_motif_info <- function(dir_motif, ncores, verbose = TRUE) {
  parfun = if (verbose) pbmcapply::pbmcmapply else parallel::mcmapply
  motif_ids = unlist(strsplit(list.files(dir_motif), '\\.cb'))
  num_inner_motifs = parfun(function(motif_id) {
    lines = readLines(paste0(motif_id, '.cb'))
    sum(grep("^>", lines))
  }, motif_ids, mc.cores = ncores)
  meta = dplyr::tibble(motif = motif_ids, num_inner_motif = num_inner_motifs)
  return(meta)
}
