#' Calculate TF-Gene P-value Matrix
#'
#' This function calculates the P-value matrix for transcription factors (TFs) and genes indicating the regulatory relationship.
#'
#' @param feather_file Path to the feather file containing Cluster-Buster scores.
#' @param fasta_consenesus_peaks Path to the FASTA file of ATAC-seq consensus peaks.
#' @param dir_null_cbscore Directory containing null Cluster-Buster scores.
#' @param peak2gene A matrix indicating the connection between ATAC-seq consensus peaks and genes.
#' @param same_chr Logical indicating whether to select control regions from the same chromosome (TRUE) or any chromosome (FALSE).
#' @param num_control_per_peak Number of matched control regions per peak. Default is 1000.
#' @param rules A vector of rules for selecting motifs. Default is `c(0.05, 0.01, 0.25)`.
#' @param dir_save Directory to save intermediate outputs, which will be loaded when the task is resumed.
#' @param resume Logical indicating whether to resume the unfinished task.
#' @param ncores Number of cores for parallel computation.
#'
#' @return A list containing (1) a P-value matrix for all TFs to all genes, where 2 indicates TF has no motifs or gene has no peaks,
#' (2) a compact P-value matrix removing the TFs and genes with P-values of 2, (3) a compact Q-value matrix where Q-values are calculated for each TF.
#' @export
#'
TF2gene <- function(feather_file, fasta_consenesus_peaks, dir_null_cbscore, peak2gene, dir_save, resume, ncores, same_chr, num_control_per_peak = 1000, rules = c(0.05, 0.01, 0.25)) {
  # Read (motifs, peaks) Cluster-Buster score matrix
  cbscore = read_cbscore(feather_file)

  # Align consensus peaks
  stopifnot(setequal(colnames(cbscore), rownames(peak2gene)))
  peak2gene = peak2gene[colnames(cbscore), ]

  # Exclude dimer and TF-pair motifs
  dimers = grep('dimer', rownames(cbscore), value = T)
  tf_pairs = grep('pair', rownames(cbscore), value = T)
  cbscore = cbscore[setdiff(rownames(cbscore), c(dimers, tf_pairs)), ]

  # Extract matched control regions for each ATAC-seq consensus peak, based on matched residue frequencies
  # mapping_mat: (consensus peaks, control regions) sparse matrix
  if (resume) {
    mapping_mat = readRDS(file.path(dir_save, 'mapping_mat.rds'))
  } else {
    meta_consensus_peak = summarize_features_of_regions(fasta = fasta_consenesus_peaks, ncores = ncores, meta_gene = TF2gene::grch38)
    mapping_mat = extract_matched_control_regions(
      meta_consensus_peak = meta_consensus_peak,
      meta_control_region = TF2gene::meta_control_region,
      n = num_control_per_peak,
      same_chr = same_chr,
      ncores = ncores
    )
    saveRDS(mapping_mat, file = file.path(dir_save, 'mapping_mat.rds'))
  }

  # Transform raw Cluster-Buster score matrix into P-value matrix based on empirical null from matched control regions
  if (resume) {
    pv_cbscore = readRDS(file.path(dir_save, 'pv_cbscore.rds'))
  } else {
    pv_cbscore = calc_pv_cbscore(
      cbscore = cbscore,
      mapping_mat = mapping_mat,
      dir_null_cbscore = dir_null_cbscore,
      n_control = num_control_per_peak,
      ncores = ncores
    )
    saveRDS(pv_cbscore, file = file.path(dir_save, 'pv_cbscore.rds'))
  }

  # Select motifs based empirical P-value matrix
  if (resume) {
    selected_motifs = readRDS(file.path(dir_save, 'selected_motifs.rds'))
  } else {
    selected_motifs = select_motifs(
      pv_cbscore = pv_cbscore,
      num_control_per_peak = num_control_per_peak,
      rules = rules,
      ncores = ncores
    )
    saveRDS(selected_motifs, file = file.path(dir_save, 'selected_motifs.rds'))
  }
  selected_pv_cbscore = pv_cbscore[selected_motifs, ]

  # Select one representative motif per TF
  tf2motif = select_representative_motif_per_TF(selected_pv_cbscore = selected_pv_cbscore, mat_tf2motif = TF2gene::mat_tf2motif)

  # Extract TF-peak P-value matrix
  pv_tf2peak = as.matrix(tf2motif %*% pmax(selected_pv_cbscore, 0.5 / num_control_per_peak))
  tfs_without_motif = rownames(tf2motif)[rowSums(tf2motif) == 0]
  pv_tf2peak[tfs_without_motif, ] = 2

  # Calculate TF-gene P-value matrix by combining 3 matrices: tf2motif, selected_pv_cbscore, peak2gene
  fisher_combination_test <- function(pvals) {
    stopifnot(all(pvals > 0 & pvals <= 1))
    pchisq(-2 * sum(log(pvals)), df = 2 * length(pvals), lower.tail = FALSE)
  }

  if (resume) {
    pv_tf2gene = readRDS(file.path(dir_save, 'pv_tf2gene.rds'))
  } else {
    pv_tf2gene = pbmcapply::pbmcmapply(function(gene) {
      peaks = names(which(peak2gene[, gene] > 0))
      if (length(peaks) == 0) return(rep(2, nrow(pv_tf2peak)))
      if (length(peaks) == 1) return(pv_tf2peak[, peaks])
      mat_pval = pv_tf2peak[, peaks]
      apply(mat_pval, 1, function(pvals) {
        if (all(pvals == 2)) return(2)
        stopifnot(all(pvals <= 1))
        fisher_combination_test(pvals)
      })
    }, colnames(peak2gene), mc.cores = ncores)
    dimnames(pv_tf2gene) = list(rownames(pv_tf2peak), colnames(peak2gene))
    saveRDS(pv_tf2gene, file = file.path(dir_save, 'pv_tf2gene.rds'))
  }

  # Get a compact version without TFs that have no motifs and genes that have no peaks
  compact_tfs = rownames(tf2motif)[Matrix::rowSums(tf2motif) > 0]
  compact_genes = colnames(peak2gene)[Matrix::colSums(peak2gene) > 0]
  stopifnot(all(pv_tf2gene[setdiff(rownames(pv_tf2gene), compact_tfs), ] == 2))
  stopifnot(all(pv_tf2gene[, setdiff(colnames(pv_tf2gene), compact_genes)] == 2))
  compact_pv_tf2gene = pv_tf2gene[compact_tfs, compact_genes]
  stopifnot(all(compact_pv_tf2gene) <= 1)

  # Calculate Q-values for each TF
  compact_qv_tf2gene = t(pbmcapply::pbmcmapply(function(i) {
    pv = compact_pv_tf2gene[i, ]
    pi0 = min(2 * mean(pv >= 0.5 & pv < 1) + mean(pv == 1), 1)
    qvalue::qvalue(pv, pi0 = pi0)$qvalues
  }, seq(nrow(compact_pv_tf2gene)), mc.cores = 36))
  dimnames(compact_qv_tf2gene) = dimnames(compact_pv_tf2gene)

  return(list(pv_tf2gene = pv_tf2gene, compact_pv_tf2gene = compact_pv_tf2gene, compact_qv_tf2gene = compact_qv_tf2gene))
}
