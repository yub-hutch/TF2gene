#' Summarize Features of Genomic Regions
#'
#' This function reads a FASTA file containing sequences of genomic regions and a metadata table of gene information.
#' It calculates the distance between the centers of the regions and the nearest TSS, and the A/C/G/T frequencies of each region.
#'
#' @param fasta A character string specifying the path to the FASTA file containing sequences of consensus peaks.
#' @param ncores An integer specifying the number of cores to use for parallel processing.
#' @param meta_gene A data frame containing metadata of genes, including chromosome and TSS information (default is TF2gene::grch38).
#'
#' @return A tibble with the following columns:
#' \describe{
#'   \item{chr}{Chromosome identifier.}
#'   \item{start}{Region start position (0-based, included).}
#'   \item{end}{Region end position (0-based, excluded).}
#'   \item{dist}{Distance from region center to nearest TSS.}
#'   \item{freq_A}{Frequency of A.}
#'   \item{freq_C}{Frequency of C.}
#'   \item{freq_G}{Frequency of G.}
#'   \item{freq_T}{Frequency of T.}
#' }
#' @export
#'
summarize_features_of_regions <- function(fasta, ncores, meta_gene = TF2gene::grch38) {
  # Read in sequences of regions
  message('Loading region sequences ...')
  fa = seqinr::read.fasta(fasta, seqtype = 'DNA', as.string = F, forceDNAtolower = T)

  # Get genomic locations of consensus peaks
  message('Extracting genomic locations ...')
  regions = names(fa)
  meta = dplyr::as_tibble(t(sapply(regions, function(region) {
    tmp = strsplit(region, '[:-]')[[1]]
    c(chr = tmp[1], start = tmp[2], end = tmp[3])
  })))
  meta = meta %>% dplyr::mutate(
    start = as.integer(start),
    end = as.integer(end)
  )

  # Calculate the distance between region centers & nearest TSS
  message('Calculating distances between region centers & nearest TSS ...')
  dist = NULL
  for (chr in unique(meta_gene$chr)) {
    message(paste0('Working on chr ', chr, ' ...'))
    meta_gene_chr = meta_gene[meta_gene$chr == chr, ]
    meta_chr = meta[meta$chr == paste0('chr', chr), ]

    meta_chr$dist = parallel::mcmapply(function(center) {
      min(abs(center - meta_gene_chr$tss))
    }, (meta_chr$start + meta_chr$end) / 2, mc.cores = ncores)

    dist = rbind(dist, meta_chr)
  }

  # A/C/G/T frequencies
  message('Calculating A/C/G/T frequencies ...')
  dist_regions = paste0(dist$chr, ':', dist$start, '-', dist$end)
  freq = do.call(rbind, parallel::mcmapply(function(seq) {
    sapply(c('a', 'c', 'g', 't'), function(letter) mean(seq == letter))
  }, fa[dist_regions], SIMPLIFY = F, mc.cores = ncores))
  colnames(freq) = paste0('freq_', c('A', 'C', 'G', 'T'))
  dist = dplyr::bind_cols(dist, dplyr::as_tibble(freq))

  return(dist)
}


#' Extract Matched control regions
#'
#' This function samples control regions that match the A/C/G/T frequencies of consensus peaks.
#'
#' @param meta_consensus_peak A tibble containing metadata of consensus peaks. See \code{\link{summarize_features_of_regions}}.
#' @param meta_control_region A tibble containing metadata for control regions.
#' @param n Number of matched control regions to sample for each consensus peak.
#' @param same_chr Logical indicating whether to select control regions from the same chromosome (TRUE) or any chromosome (FALSE).
#' @param ncores An integer specifying the number of cores to use for parallel processing.
#'
#' @return A sparse matrix recording the distances between each consensus peak and its matched control regions.
#' Distance of 0 is replaced by 1e-10 for sparse matrix representation.
#' @export
extract_matched_control_regions <- function(meta_consensus_peak, meta_control_region, n, same_chr, ncores) {
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
  if (same_chr) {
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
          control_id = paste0(curr_meta_control_region$chr[top_control_indices], ":",
                              curr_meta_control_region$start[top_control_indices], "-",
                              curr_meta_control_region$end[top_control_indices]),
          distance = distances[top_control_indices]
        ))
      }
      record
    }, unique(meta_consensus_peak$chr), SIMPLIFY = FALSE, mc.cores = ncores))
  } else {
    result = do.call(rbind, pbmcapply::pbmcmapply(function(i) {
      consensus_peak = meta_consensus_peak[i, ]
      distances = sqrt((meta_control_region$freq_A - consensus_peak$freq_A)^2 +
                         (meta_control_region$freq_C - consensus_peak$freq_C)^2 +
                         (meta_control_region$freq_G - consensus_peak$freq_G)^2 +
                         (meta_control_region$freq_T - consensus_peak$freq_T)^2)
      top_control_indices = order(distances)[1:n]
      data.frame(
        consensus_id = paste0(consensus_peak$chr, ":", consensus_peak$start, "-", consensus_peak$end),
        control_id = paste0(meta_control_region$chr[top_control_indices], ":",
                            meta_control_region$start[top_control_indices], "-",
                            meta_control_region$end[top_control_indices]),
        distance = distances[top_control_indices]
      )
    }, seq(nrow(meta_consensus_peak)), SIMPLIFY = FALSE, mc.cores = ncores))
  }
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
