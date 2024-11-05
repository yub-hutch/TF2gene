#' Select Active Motifs
#'
#' Calculates enrichment proportions for each motif based on transformed Cluster-Buster scores.
#'
#' @param trans_cbscore Matrix or data frame of transformed Cluster-Buster scores (See \code{\link{transform_cbscore}}).
#'
#' @return Tibble with `motif` and `prop_enrich` columns.
#' @export
select_active_motifs <- function(trans_cbscore) {
  thrs = seq(0.5, 0.9, by = 0.1)
  props_enrich = sapply(thrs, function(thr) {
    prop_above_thr = rowMeans(trans_cbscore > thr)
    prop_below_thr_under_null = thr / (1 - thr) * prop_above_thr
    1 - prop_above_thr - prop_below_thr_under_null
  })
  min_prop_enrich = apply(props_enrich, 1, min)
  dplyr::tibble(motif = names(min_prop_enrich), min_prop_enrich = min_prop_enrich)
}
