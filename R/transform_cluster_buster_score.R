#' Transform Cluster-Buster Scores to P-values
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
  # Check arguments
  stopifnot(setequal(rownames(mapping_mat), colnames(cbscore)))
  consensus_peaks = rownames(mapping_mat)

  # Subset motifs
  message('Subsetting motifs ...')
  cbscore = cbscore[c(depleted_motifs, enriched_motifs), consensus_peaks]

  # Set score of depleted motifs to -1
  cbscore[depleted_motifs, ] = -1

  # Extract top matched control regions
  message('Extract top matched control regions ...')
  control_regions = extract_top_matched_control_regions(mapping_mat, n_control, ncores)
  stopifnot(identical(names(control_regions), consensus_peaks))
  rm(mapping_mat, n_control)
  gc()

  # Transform Cluster-Buster score of enriched motifs
  message('Transforming Cluster-Buster score of enriched motifs ...')
  vs_cbscore = lapply(enriched_motifs, function(motif) cbscore[motif, ])
  cbscore[enriched_motifs, ] = t(pbmcapply::pbmcmapply(function(motif, v_cbscore) {
    all_null_score = read_null_cbscore(motif, dir_null_cbscore)
    pvs = mapply(function(regions, score) {
      matched_null_score = all_null_score[regions]
      mean(matched_null_score >= score)
    }, control_regions, v_cbscore)
    pvs
  }, enriched_motifs, vs_cbscore, mc.cores = ncores))

  return(cbscore)
}
