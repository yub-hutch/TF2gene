#' Load Null Cluster-Buster Scores
#'
#' This function loads the null cluster-buster scores for a given motif.
#'
#' @param motif A character string representing the motif name.
#' @param dir_null_cbscore A character string representing the directory containing the null cluster-buster scores.
#' @return A named numeric vector of scores.
#' @export
load_null_cbscore <- function(motif, dir_null_cbscore) {
  fname = file.path(dir_null_cbscore, paste0(motif, '.motifs_vs_regions.scores.feather'))
  raw = arrow::read_feather(fname)
  score = setNames(raw[[1]], raw$regions)
  return(score)
}


#' Load Cluster-Buster Scores
#'
#' This function loads the cluster-buster scores from a feather file.
#'
#' @param feather A character string representing the path to the feather file.
#' @return A matrix of scores with motifs as rows and regions as columns.
#' @export
load_cbscore <- function(feather) {
  raw = arrow::read_feather(feather)
  mat = t(as.matrix(raw[, setdiff(names(raw), 'regions')]))
  colnames(mat) = raw$regions
  mat = mat[, grepl("^chr([1-9]|1[0-9]|2[0-2]|X|Y):", colnames(mat))]
  return(mat)
}


#' Select Motifs Based on Cluster-Buster Scores
#'
#' This function selects motifs based on their cluster-buster scores and control peaks.
#'
#' @param cbscore A matrix of cluster-buster scores.
#' @param control_peaks A character vector of control peak names.
#' @param dir_null_cbscore A character string representing the directory containing the null cluster-buster scores.
#' @param ncores An integer specifying the number of cores to use for parallel processing.
#' @return A tibble with motifs, log fold changes, and p-values.
#' @export
select_motifs <- function(cbscore, control_peaks, dir_null_cbscore, ncores) {
  cbscore_list = sapply(rownames(cbscore), simplify = F, function(motif) cbscore[motif, ])
  do.call(rbind, pbmcapply::pbmcmapply(function(motif, score) {
    # Set matched null cluster-buster score
    all_null_score = load_null_cbscore(motif, dir_null_cbscore)
    matched_null_score = all_null_score[control_peaks]

    # Fold change & P-value
    logfc = log(mean(score) / mean(matched_null_score))
    pv = wilcox.test(score, matched_null_score, alternative = 'greater')$p.value

    dplyr::tibble(motif = motif, logfc = logfc, pv = pv)
  }, names(cbscore_list), cbscore_list, SIMPLIFY = F, mc.cores = ncores))
}
