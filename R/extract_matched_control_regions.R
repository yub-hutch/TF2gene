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



#' Extract Best Matched Control Regions
#'
#' This function finds the best-matched control region for each consensus peak.
#'
#' @param mapping_mat A sparse matrix where rows represent consensus peaks and columns
#'   represent control regions. Each entry in the matrix is a value indicating the
#'   distance between the consensus peak and the control region.
#'   See \code{\link{extract_matched_control_regions}}.
#' @param ncores An integer specifying the number of cores to use for parallel processing.
#'
#' @return A named vector of the best control region for each consensus peak.
#'
#' @export
extract_best_matched_control_region <- function(mapping_mat, ncores) {
  control_regions = colnames(mapping_mat)
  chunk_indices = split(seq(nrow(mapping_mat)), cut(seq(nrow(mapping_mat)), ncores))
  mat_chunks = lapply(chunk_indices, function(indices) mapping_mat[indices, , drop = FALSE])
  result = do.call(rbind, pbmcapply::pbmcmapply(function(sub_mat) {
    curr = NULL
    for (i in seq(nrow(sub_mat))) {
      v = sub_mat[i, ]
      non_zero_indices = which(v != 0)
      min_nonzero_index = non_zero_indices[which.min(v[non_zero_indices])]
      best_control_region = data.frame(control_region = control_regions[min_nonzero_index], distance = min(v[non_zero_indices]))
      curr = rbind(curr, best_control_region)
    }
    curr
  }, mat_chunks, mc.cores = ncores, SIMPLIFY = F))
  result$consensus_peak = rownames(mapping_mat)
  result = dplyr::as_tibble(result)[, c('consensus_peak', 'control_region', 'distance')]
  return(result)
}


#' Extract Top Matched Control Regions
#'
#' This function finds the top-matched control regions for each consensus peak.
#'
#' @param mapping_mat A sparse matrix where rows represent consensus peaks and columns
#'   represent control regions. Each entry in the matrix is a value indicating the
#'   distance between the consensus peak and the control region.
#'   See \code{\link{extract_matched_control_regions}}.
#' @param n_control Number of matched control regions to extract for each consensus peak.
#' @param ncores An integer specifying the number of cores to use for parallel processing.
#'
#' @return A named list of the top matched control regions for each consensus peak.
#'
#' @export
extract_top_matched_control_regions <- function(mapping_mat, n_control, ncores) {
  chunk_indices = split(seq(nrow(mapping_mat)), cut(seq(nrow(mapping_mat)), ncores))
  mat_chunks = lapply(chunk_indices, function(indices) mapping_mat[indices, , drop = FALSE])

  Reduce(c, pbmcapply::pbmcmapply(function(sub_mat) {
    curr = list()
    consensus_peaks = rownames(sub_mat)
    for (consensus_peak in consensus_peaks) {
      mapping = sub_mat[consensus_peak, ]
      mapping = mapping[mapping > 0]
      stopifnot(length(mapping) >= n_control)
      control_regions = names(head(sort(mapping), n_control))
      curr[[consensus_peak]] = control_regions
    }
    curr
  }, mat_chunks, mc.cores = ncores, SIMPLIFY = F))
}
