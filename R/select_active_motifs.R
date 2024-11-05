#' Cauchy Combination Test
#'
#' This function combines \( p \)-values using the Cauchy Combination Test, which is effective
#' for handling dependencies among tests.
#'
#' @param pvals Numeric vector of p-values (each between 1e-16 and 1).
#' @return Combined p-value from the Cauchy test.
#' @export
cauchy_combination_test <- function(pvals) {
  stopifnot(all(!is.na(pvals)) & all(pvals > 1e-16) & all(pvals < 1))
  weights = rep(1 / length(pvals), length(pvals))
  stat = sum(weights * tan((0.5 - pvals) * pi))
  ifelse(stat > 1e15, 1 / stat / pi, 1 - pcauchy(stat))
}


#' Select Active Motifs Based on Transformed Cluster-Buster Scores
#'
#' This function applies the Cauchy Combination Test to transformed Cluster-Buster scores and
#' returns a tibble of motifs and combined p-values.
#'
#' @param trans_cbscore Matrix of transformed Cluster-Buster scores (rows are motifs).
#' @param n_control Numeric, the number of control regions used.
#' @return A tibble with motifs and their combined p-values.
#' @export
select_active_motifs <- function(trans_cbscore, n_control) {
  epsilon = 0.1 / n_control
  trans_cbscore[trans_cbscore == 0] = epsilon
  trans_cbscore[trans_cbscore == 1] = 1 - epsilon
  pvals = apply(trans_cbscore, 1, cauchy_combination_test)
  dplyr::tibble(motif = rownames(trans_cbscore), pv = pvals)
}
