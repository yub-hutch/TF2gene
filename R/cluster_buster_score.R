#' Read Cluster-Buster Scores
#'
#' This function reads the cluster-buster scores from a feather file.
#'
#' @param feather A character string representing the path to the feather file.
#'
#' @return A matrix of scores with motifs as rows and regions as columns.
#' @export
#'
read_cbscore <- function(feather) {
  raw = arrow::read_feather(feather)
  mat = t(as.matrix(raw[, setdiff(names(raw), 'regions')]))
  colnames(mat) = raw$regions
  mat = mat[, grepl("^chr([1-9]|1[0-9]|2[0-2]|X|Y):", colnames(mat))]
  return(mat)
}


#' Read Null Cluster-Buster Scores
#'
#' This function reads the cluster-buster scores on control regions for a given motif.
#'
#' @param motif A character string representing the motif name.
#' @param dir_null_cbscore A character string representing the directory containing the null cluster-buster scores.
#'
#' @return A named numeric vector of scores.
#' @export
#'
read_null_cbscore <- function(motif, dir_null_cbscore) {
  fname = file.path(dir_null_cbscore, paste0(motif, '.motifs_vs_regions.scores.feather'))
  raw = arrow::read_feather(fname)
  score = setNames(raw[[1]], raw$regions)
  return(score)
}


#' Transform Cluster-Buster Scores to P-values
#'
#' Transforms Cluster-Buster scores into P-values by comparing to (motif, consensus peak)-specific null scores.
#'
#' @param cbscore Matrix of Cluster-Buster scores with motifs as rows and consensus peaks as columns.
#' @param mapping_mat Sparse matrix mapping consensus peaks to control regions. See \code{\link{extract_matched_control_regions}}.
#' @param dir_null_cbscore Directory path with pre-computed null distributions for each motif.
#' @param n_control Number of control regions used in null score calculation in alignment with mapping_mat.
#' @param ncores Number of cores for parallel processing.
#'
#' @return A P-value matrix.
#' @export
#'
calc_pv_cbscore <- function(cbscore, mapping_mat, dir_null_cbscore, n_control, ncores) {
  # Check arguments & format inputs
  message('Formatting inputs ...')
  consensus_peaks = rownames(mapping_mat)
  control_regions = colnames(mapping_mat)
  mapping_mat = mapping_mat > 0

  stopifnot(setequal(consensus_peaks, colnames(cbscore)))
  cbscore = cbscore[, consensus_peaks]
  motifs = rownames(cbscore)

  detected_n_controls = Matrix::rowSums(mapping_mat)
  stopifnot(all(detected_n_controls == n_control))

  # Transform Cluster-Buster score
  message('Splitting cbscore for parallel computation ...')
  vs_cbscore = lapply(motifs, function(motif) cbscore[motif, ])

  message('Transforming Cluster-Buster score ...')
  pv = t(pbmcapply::pbmcmapply(function(motif, v_cbscore) {
    all_null_score = read_null_cbscore(motif, dir_null_cbscore)[control_regions]

    # To do efficient sparse matrix operation, otherwise computationally infeasible
    diag_all_null_score = Matrix::Diagonal(x = all_null_score + 1) # Add 1 to distinguish blank elements and 0 score
    mapping_mat_with_null_cbscore = mapping_mat %*% diag_all_null_score
    comparison = Matrix::sparseMatrix(
      i = mapping_mat_with_null_cbscore@i + 1, # @i is 0-based
      j = Matrix::summary(mapping_mat_with_null_cbscore)$j, # Matrix::summary()$j is 1-based
      x = mapping_mat_with_null_cbscore@x >= (v_cbscore[mapping_mat_with_null_cbscore@i + 1] + 1),
      dims = dim(mapping_mat_with_null_cbscore)
    )
    Matrix::rowSums(comparison) / n_control
  }, motifs, vs_cbscore, mc.cores = ncores))
  dimnames(pv) = dimnames(cbscore)

  return(pv)
}

