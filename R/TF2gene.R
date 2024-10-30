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
  keep = (peak2gene_distance@x > 0) & (peak2gene_distance@x < distance_thr)
  peak2gene = Matrix::sparseMatrix(
    i = peak2gene_distance@i[keep] + 1, # @i is 0-based
    j = Matrix::summary(peak2gene_distance)$j[keep], # summary()$j is 1-based
    x = 1,
    dims = dim(peak2gene_distance),
    dimnames = dimnames(peak2gene_distance)
  )

  # Calculate TF-gene score
  tf2gene = as.matrix(tf2peak %*% peak2gene)
  return(tf2gene)
}
