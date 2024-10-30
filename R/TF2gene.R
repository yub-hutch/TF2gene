#' Calculate TF-Gene Score
#'
#' Computes the transcription factor (TF)-to-gene score by combining TF-peak
#' scores with peak-to-gene distances, using a specified distance threshold.
#'
#' @param tf2peak A matrix of TF-to-peak scores.
#' @param peak2gene_distance A sparse matrix with peak-to-gene distances.
#' @param distance_thr Numeric. The maximum peak-to-gene distance threshold.
#' @return A matrix representing TF-to-gene scores.
#' @export
#'
calc_tf2gene_score <- function(tf2peak, peak2gene_distance, distance_thr) {
  # Check arguments
  consensus_peaks = rownames(peak2gene_distance)
  stopifnot(setequal(colnames(tf2peak), consensus_peaks))
  tf2peak = tf2peak[, consensus_peaks]

  # Subset connections between consensus peak and genes
  peak2gene = (peak2gene_distance > 0) & (peak2gene_distance < peak2gene_distance_thr)

  # Calculate TF-gene score
  tf2gene = tf2peak %*% peak2gene
  if (Matrix::mean(tf2gene == 0) < 0.2) tf2gene = as.matrix(tf2gene)
  return(tf2gene)
}
