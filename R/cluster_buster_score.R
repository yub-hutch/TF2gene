#' Read Null Cluster-Buster Scores
#'
#' This function reads the cluster-buster scores on control regions for a given motif.
#'
#' @param motif A character string representing the motif name.
#' @param dir_null_cbscore A character string representing the directory containing the null cluster-buster scores.
#' @return A named numeric vector of scores.
#' @export
read_null_cbscore <- function(motif, dir_null_cbscore) {
  fname = file.path(dir_null_cbscore, paste0(motif, '.motifs_vs_regions.scores.feather'))
  raw = arrow::read_feather(fname)
  score = setNames(raw[[1]], raw$regions)
  return(score)
}


#' Read Cluster-Buster Scores
#'
#' This function reads the cluster-buster scores from a feather file.
#'
#' @param feather A character string representing the path to the feather file.
#' @return A matrix of scores with motifs as rows and regions as columns.
#' @export
read_cbscore <- function(feather) {
  raw = arrow::read_feather(feather)
  mat = t(as.matrix(raw[, setdiff(names(raw), 'regions')]))
  colnames(mat) = raw$regions
  mat = mat[, grepl("^chr([1-9]|1[0-9]|2[0-2]|X|Y):", colnames(mat))]
  return(mat)
}
