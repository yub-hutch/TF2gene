#' Select Active Motifs
#'
#' Calculates enrichment proportions for each motif based on transformed Cluster-Buster scores.
#'
#' @param trans_cbscore Matrix or data frame of transformed Cluster-Buster scores (See \code{\link{transform_cbscore}}).
#'
#' @return Tibble with `motif` and `prop_enrich` columns.
#' @export
select_active_motifs <- function(trans_cbscore) {
  prop_enrich = 1 - 2 * rowMeans(trans_cbscore > 0.5)
  dplyr::tibble(motif = names(prop_enrich), prop_enrich = prop_enrich)
}
