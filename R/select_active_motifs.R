#' Select Active Motifs Based on Cluster-Buster Scores
#'
#' This function selects active motifs by comparing their cluster-buster scores on consensus peaks and matched control regions.
#'
#' @param cbscore A matrix of cluster-buster scores of consensus peaks.
#' @param control_regions A names vector of best control regions for each consensus peak.
#' See \code{\link{extract_best_matched_control_region}}.
#' @param dir_null_cbscore A character string representing the directory containing the cluster-buster scores on control regions.
#' @param ncores An integer specifying the number of cores to use for parallel processing.
#' @param dir_out A character string representing the directory to save the results.
#' @return This function is called by its side effect.
#' @export
select_active_motifs_with_cbscore <- function(cbscore, control_regions, dir_null_cbscore, ncores, dir_out) {
  if (!dir.exists(dir_out)) dir.create(dir_out)

  if (!(setequal(names(control_regions), colnames(cbscore)) &
        (ncol(cbscore) == length(control_regions)))) stop("cbscore & control_regions don't match")

  # Run
  message('Preparing to run ...')
  cbscore_list = sapply(rownames(cbscore), simplify = F, function(motif) cbscore[motif, ])
  rm(cbscore)
  gc()

  unfinished_motifs = names(cbscore_list)
  message(paste0(length(unfinished_motifs), ' motifs detected'))

  # Compare
  message('Comparing Cluster-Buster scores on consensus peaks & matched control regions ...')
  while (length(unfinished_motifs) > 0) {
    junk = tryCatch({
      pbmcapply::pbmcmapply(function(motif, score) {
        # Set matched null cluster-buster score
        all_null_score = read_null_cbscore(motif, dir_null_cbscore)
        matched_null_score = all_null_score[control_regions]

        # Fold change, P-value, & AUC
        logfc = log(mean(score) / mean(matched_null_score))
        pv = wilcox.test(score, matched_null_score)$p.value
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

