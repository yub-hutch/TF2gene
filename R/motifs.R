#' Select Motifs Based on P-value and Q-value Rules
#'
#' This function selects motifs based on empirical P-value and Q-value rules.
#'
#' @param pv_cbscore A matrix of P-values for motifs to ATAC-seq consensus peaks.
#' @param num_control_per_peak Number of control regions per peak.
#' @param rules A vector of rules for selecting motifs.
#' @param ncores Number of cores to use for parallel processing.
#'
#' @return A vector of selected motif names.
#' @export
#'
select_motifs <- function (pv_cbscore, num_control_per_peak, rules, ncores) {
  qv_cbscore = t(pbmcapply::pbmcmapply(function(i) {
    pv = pmax(pv_cbscore[i, ], 0.5 / num_control_per_peak) # minimum nonzero empirical P-value is 1/num_control_per_peak
    pi0 = min(2 * mean(pv >= 0.5 & pv < 1) + mean(pv == 1), 1)
    qvalue::qvalue(pv, pi0 = pi0)$qvalues
  }, seq(nrow(pv_cbscore)), mc.cores = ncores))
  dimnames(qv_cbscore) = dimnames(pv_cbscore)

  rule1 = rowMeans(pv_cbscore < rules[1]) > rules[1]
  rule2 = rowMeans(pv_cbscore < rules[2]) > rules[2]
  rule3 = apply(qv_cbscore, 1, min) < rules[3]

  selected_motifs = rownames(pv_cbscore)[rule1 & rule2 & rule3]
  return(selected_motifs)
}


#' Select Representative Motif per Transcription Factor
#'
#' This function selects one representative motif per transcription factor (TF) based on various criteria.
#'
#' @param selected_pv_cbscore A P-value matrix for selected motifs to consensus peaks.
#' @param mat_tf2motif A matrix indicating the connection between TFs and motifs. Default is `TF2gene::mat_tf2motif`.
#' @return A matrix indicating the selected representative motif for each TF.
#' @export
#'
select_representative_motif_per_TF <- function(selected_pv_cbscore, mat_tf2motif = TF2gene::mat_tf2motif) {
  # Only keep connections with direct annotation
  tf2motif = TF2gene::mat_tf2motif >= 3
  tf2motif = tf2motif[, rownames(selected_pv_cbscore)]

  message(message('# TFs with connection: ', sum(Matrix::rowSums(tf2motif) > 0)))
  message(message('# Connections: ', sum(tf2motif@x)))

  # Select one representative motif per TF
  for(tf in rownames(tf2motif)) {
    motifs = names(which(tf2motif[tf, ]))
    if (length(motifs) < 2) next
    if (length(motifs) >= 3) {
      # For TFs with >= 3 motifs, select the motif with the strongest correlations with the other motifs
      cor_mat = cor(t(selected_pv_cbscore[motifs, ]), method = 'spearman')
      diag(cor_mat) = NA
      median_cors = apply(cor_mat, 1, median, na.rm = T)
      selected_motif = motifs[which.max(median_cors)]
    }
    if (length(motifs) == 2) {
      # For TFs with only 2 motifs,
      # 1. If one is singlet, one is metacluster, choose the metacluster
      # 2. If both are metaclusters, choose the one with more member motifs
      # 3. If both are singlets, or metaclusters with same number of member motifs, choose the one with P-value < 0.05 on more consensus peaks
      is_meta = substr(motifs, 1, 11) == 'metacluster'
      if (sum(is_meta) == 1) selected_motif = motifs[is_meta]
      if (sum(is_meta) == 2) {
        n_member = TF2gene::meta_motif$num_inner_motif[match(motifs, TF2gene::meta_motif$motif)]
        if (n_member[1] != n_member[2]) {
          selected_motif = motifs[which.max(n_member)]
        } else {
          prop_sig = rowMeans(selected_pv_cbscore[motifs, ] < 0.05)
          selected_motif = motifs[which.max(prop_sig)]
        }
      }
      if (sum(is_meta) == 0) {
        prop_sig = rowMeans(selected_pv_cbscore[motifs, ] < 0.05)
        selected_motif = motifs[which.max(prop_sig)]
      }
    }
    tf2motif[tf, setdiff(colnames(tf2motif), selected_motif)] = F
  }
  return(tf2motif)
}
