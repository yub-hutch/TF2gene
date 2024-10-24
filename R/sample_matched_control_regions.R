#' Extract Matched control regions
#'
#' This function samples control regions that match the A/C/G/T frequencies of consensus peaks.
#'
#' @param meta_consensus_peak A tibble containing metadata of consensus peaks. See \code{\link{summarize_features_of_regions}}.
#' @param meta_control_region A tibble containing metadata for control regions.
#' @param n Number of matched control regions to sample for each consensus peak.
#' @param ncores An integer specifying the number of cores to use for parallel processing.
#'
#' @return A sparse matrix recording the distances between each consensus peak and its matched control regions.
#' @export
extract_matched_control_regions <- function(meta_consensus_peak, meta_control_region, n, ncores) {
  message(paste0(nrow(meta_consensus_peak), ' consensus peaks detected.'))
  message(paste0(nrow(meta_control_region), ' controls regions detected.'))

  # Exclude control regions that overlap with consensus peaks
  message('Excluding control regions that overlap with consensus peaks ...')
  meta_control_region = do.call(rbind, parallel::mcmapply(function(chr) {
    curr_meta_consensus_peak = meta_consensus_peak[meta_consensus_peak$chr == chr, ]
    curr_meta_control_region = meta_control_region[meta_control_region$chr == chr, ]
    ir_meta_consensus_peak = IRanges::IRanges(
      start = curr_meta_consensus_peak$start,
      end = curr_meta_consensus_peak$end
    )
    ir_meta_control_region = IRanges::IRanges(
      start = curr_meta_control_region$start,
      end = curr_meta_control_region$end
    )
    overlaps = IRanges::findOverlaps(ir_meta_control_region, ir_meta_consensus_peak)
    if (length(overlaps) > 0) curr_meta_control_region[-S4Vectors::queryHits(overlaps), ] else curr_meta_control_region
  }, unique(meta_consensus_peak$chr), SIMPLIFY = F, mc.cores = ncores))
  message(paste0(nrow(meta_control_region), ' controls peaks remaining after excluding overlaps with consensus peaks.'))

  # Extract matched control regions
  message('Extracting matched control regions ...')
  result = do.call(rbind, pbmcapply::pbmcmapply(function(chr) {
    curr_meta_consensus_peak = meta_consensus_peak[meta_consensus_peak$chr == chr, ]
    curr_meta_control_region = meta_control_region[meta_control_region$chr == chr, ]
    record = NULL
    for (i in 1:nrow(curr_meta_consensus_peak)) {
      consensus_peak = curr_meta_consensus_peak[i, ]
      distances = sqrt((curr_meta_control_region$freq_A - consensus_peak$freq_A)^2 +
                         (curr_meta_control_region$freq_C - consensus_peak$freq_C)^2 +
                         (curr_meta_control_region$freq_G - consensus_peak$freq_G)^2 +
                         (curr_meta_control_region$freq_T - consensus_peak$freq_T)^2)
      top_control_indices = order(distances)[1:n]
      record = rbind(record, data.frame(
        consensus_id = paste0(consensus_peak$chr, ":", consensus_peak$start, "-", consensus_peak$end),
        control_id = paste0(curr_meta_control_region$chr[top_control_indices], ":", curr_meta_control_region$start[top_control_indices],
                            "-", curr_meta_control_region$end[top_control_indices]),
        distance = distances[top_control_indices]
      ))
    }
    record
  }, unique(meta_consensus_peak$chr), SIMPLIFY = FALSE, mc.cores = ncores))
  result$distance[result$distance == 0] = 1e-10

  # Construct the sparse matrix
  message('Buiding sparse matrix ...')
  consensus_peaks = paste0(meta_consensus_peak$chr, ":", meta_consensus_peak$start, "-", meta_consensus_peak$end)
  control_regions = unique(result$control_id)
  mat = Matrix::sparseMatrix(
    i = match(result$consensus_id, consensus_peaks),
    j = match(result$control_id, control_regions),
    x = result$distance,
    dims = c(length(consensus_peaks), length(control_regions)),
    dimnames = list(consensus_peaks, control_regions)
  )
  return(mat)
}
