#' @import Matrix
NULL


#' Transform Cluster-Buster Scores
#'
#' Transforms Cluster-Buster scores for enriched and depleted motifs.
#'
#' @param cbscore Matrix of Cluster-Buster scores with motifs as rows and consensus peaks as columns.
#' @param depleted_motifs Vector of motif names set as depleted, with scores set to zero.
#' @param enriched_motifs Vector of motif names to be enriched and transformed.
#' @param mapping_mat Sparse matrix mapping consensus peaks to control regions. See \code{\link{extract_matched_control_regions}}.
#' @param dir_null_cbscore Directory path with precomputed null distributions for each motif.
#' @param n_control Number of control regions used in null score calculation.
#' @param ncores Number of cores for parallel processing.
#'
#' @return Modified `cbscore` matrix with transformed scores.
#' @export
#'
transform_cbscore <- function(cbscore, depleted_motifs, enriched_motifs, mapping_mat, dir_null_cbscore, n_control, ncores) {
  message('Subsetting motifs ...')
  cbscore = cbscore[c(depleted_motifs, enriched_motifs), ]

  # Set score of depleted motifs to 0
  cbscore[depleted_motifs, ] = 0

  # Transform Cluster-Buster score of enriched motifs
  message('Transforming Cluster-Buster score of enriched motifs ...')
  consensus_peaks = colnames(cbscore)
  cbscore[enriched_motifs, ] = t(pbmcapply::pbmcmapply(function(motif, v_cbscore) {
    all_null_score = read_null_cbscore(motif, dir_null_cbscore)
    pvs = setNames(rep(NA, length(consensus_peaks)), consensus_peaks)
    for (consensus_peak in consensus_peaks) {
      mapping = mapping_mat[consensus_peak, ]
      mapping = mapping[mapping > 0]
      stopifnot(length(mapping) >= n_control)
      control_regions = names(head(sort(mapping), n_control))
      matched_null_score = all_null_score[control_regions]
      pvs[consensus_peak] = mean(matched_null_score >= v_cbscore[consensus_peak])
    }
    pvs
  }, enriched_motifs, lapply(enriched_motifs, function(motif) cbscore[motif, ]), mc.cores = ncores))

  return(cbscore)
}
