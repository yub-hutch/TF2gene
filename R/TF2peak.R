#' Calculate TF-Peak Score
#'
#' Computes a transcription factor (TF)-to-peak score with transformed Cluster-Buster score matrix
#' of enriched motifs and TF-motif connection matrix.
#'
#' @param trans_cbscore A matrix of transformed Cluster-Buster scores, i.e., P-values, of enriched motifs.
#' See \code{\link{transform_cbscore}} and \code{\link{select_active_motifs_with_cbscore}}.
#' @param n_control Number of control regions used to generate \code{trans_cbscore}.
#' @param mat_tf2motif A binary matrix indicating TF-motif connections. See \code{?mat_tf2motif} for default.
#' @param tf2motif_level Integer (≤3) indicating the minimum level of TF-motif connections to retain.
#' 1: by similarity; 2: by orthology; 3. direct annotation.
#' @return A matrix representing TF-to-peak scores.
#' @export
#'
calc_tf2peak_score <- function(trans_cbscore, n_control, mat_tf2motif = mat_tf2motif, tf2motif_level) {
  # Check arguments
  stopifnot(all(rownames(trans_cbscore) %in% colnames(mat_tf2motif)))
  stopifnot(tf2motif_level <= 3)
  message(paste0(nrow(mat_tf2motif), ' TFs detected'))

  # Subset TF-motif connection
  message('Subsetting TF-motif connection above the specified level ...')
  mat_tf2motif = (mat_tf2motif >= tf2motif_level)

  message('Removing TFs without connection ...')
  motifs = rownames(trans_cbscore)
  mat_tf2motif = mat_tf2motif[, motifs]
  keep_TF = Matrix::rowSums(mat_tf2motif) > 0
  stopifnot(sum(keep_TF) > 0)
  mat_tf2motif = mat_tf2motif[keep_TF, ]
  message(paste0(nrow(mat_tf2motif), ' TFs left'))

  # Transform P-value
  message('Transforming the input P-value matrix ...')
  epsilon = 0.1 / n_control
  trans_cbscore[trans_cbscore == 0] = epsilon
  motif2peak = -log10(trans_cbscore)

  # Calculate TF-peak score
  message('Calculating TF-peak score ...')
  tf2peak = mat_tf2motif %*% motif2peak
  if (Matrix::mean(tf2peak == 0) < 0.2) tf2peak = as.matrix(tf2peak)
  return(tf2peak)
}
