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
#' @param dir_out A character string representing the directory to save the results.
#' @return This function is called by its side effect.
#' @export
select_motifs_with_cbscore <- function(cbscore, control_peaks, dir_null_cbscore, ncores, dir_out) {
  if (!dir.exists(dir_out)) dir.create(dir_out)

  # Run
  message('Preparing to run ...')
  cbscore_list = sapply(rownames(cbscore), simplify = F, function(motif) cbscore[motif, ])
  rm(cbscore)
  gc()

  unfinished_motifs = names(cbscore_list)
  message(paste0(length(unfinished_motifs), ' motifs detected'))

  while (length(unfinished_motifs) > 0) {
    junk = tryCatch({
      pbmcapply::pbmcmapply(function(motif, score) {
        # Set matched null cluster-buster score
        all_null_score = load_null_cbscore(motif, dir_null_cbscore)
        matched_null_score = all_null_score[control_peaks]

        # Fold change, P-value, & AUC
        logfc = log(mean(score) / mean(matched_null_score))
        pv = wilcox.test(score, matched_null_score, alternative = 'greater')$p.value
        auc = as.numeric(pROC::auc(
          response = c(rep(1, length(score)), rep(0, length(matched_null_score))),
          predictor = c(score, matched_null_score),
          levels = c(0, 1),
          direction = '<'
        ))

        res = dplyr::tibble(motif = motif, logfc = logfc, pv = pv, auc = auc)

        # Save
        saveRDS(res, file = file.path(dir_out, paste0(motif, '.rds')))

        return(NULL)
      }, names(cbscore_list), cbscore_list, SIMPLIFY = F, mc.cores = ncores)
    }, error = function(e) 'unfinished')

    finished_motifs = sapply(strsplit(list.files(dir_out), '.rds'), head, n = 1)
    unfinished_motifs = setdiff(names(cbscore_list), finished_motifs)

    if (length(unfinished_motifs) > 0) {
      message(paste0(length(unfinished_motifs), ' motifs unfinished. Retry ...'))
      cbscore_list = cbscore_list[unfinished_motifs]
    } else {
      message('done.')
    }
  }
  return(NULL)
}


#' Collect Motif Selection Results
#'
#' This function collects outputs of \link{\code{select_motifs_with_cbscore}} from a specified directory.
#'
#' @param dir_out A character string specifying the directory containing the RDS files.
#' @param ncores An integer specifying the number of cores to use for parallel processing.
#'
#' @return A tibble with columns motif, logfc, and pv.
#' @export
collect_motif_selection <- function(dir_out, ncores) {
  do.call(rbind, parallel::mcmapply(function(file) {
    readRDS(file.path(dir_out, file))
  }, list.files(dir_out), SIMPLIFY = F, mc.cores = ncores))
}
