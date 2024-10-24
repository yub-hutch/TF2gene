#' Sample Matched control regions
#'
#' This function samples control regions that match the A/C/G/T frequencies of consensus peaks.
#'
#' @param meta_peak A tibble containing metadata for consensus peaks. See \code{\link{summarize_consensus_peaks}}.
#' @param meta_control_region A tibble containing metadata for control regions.
#' @param nbins Number of bins in each dimension. May be a scalar or a 2 element vector.
#' @param n_control A numeric value specifying the number of control regions to sample.
#' @param ncores An integer specifying the number of cores to use for parallel processing.
#' @param seed A numeric value specifying the seed for random sampling.
#'
#' @return A tibble of subset of meta_control_region containing the sampled control regions.
#' @export
sample_matched_control_regions <- function(meta_peak, meta_control_region, nbins, n_control, ncores, seed) {
  set.seed(seed)

  message(paste0(nrow(meta_peak), ' consensus peaks detected.'))
  message(paste0(nrow(meta_control_region), ' controls peaks detected.'))

  # Exclude control regions that overlap with consensus peaks
  meta_control_region = do.call(rbind, parallel::mcmapply(function(chr) {
    curr_meta_peak = meta_peak[meta_peak$chr == chr, ]
    curr_meta_control_region = meta_control_region[meta_control_region$chr == chr, ]
    ir_meta_peak = IRanges::IRanges(
      start = curr_meta_peak$start,
      end = curr_meta_peak$end
    )
    ir_meta_control_region = IRanges::IRanges(
      start = curr_meta_control_region$start,
      end = curr_meta_control_region$end
    )
    overlaps = IRanges::findOverlaps(ir_meta_control_region, ir_meta_peak)
    if (length(overlaps) > 0) curr_meta_control_region[-S4Vectors::queryHits(overlaps), ] else curr_meta_control_region
  }, unique(meta_peak$chr), SIMPLIFY = F, mc.cores = ncores))
  message(paste0(nrow(meta_control_region), ' controls peaks remaining after excluding overlaps with consensus peaks.'))

  # Exclude consensus peaks with features falling outside those ranges of consensus peaks
  meta_control_region = meta_control_region %>%
    dplyr::filter(dist >= min(meta_peak$dist),
                  dist <= max(meta_peak$dist),
                  gc_content >= min(meta_peak$gc_content),
                  gc_content <= max(meta_peak$gc_content))
  message(paste0(nrow(meta_control_region), ' controls peaks remaining after excluding out-of-range ones.'))

  # Fit density
  message('Fitting density for consensus peaks ...')
  fit = gplots::hist2d(
    x = log10(1 + meta_peak$dist),
    y = meta_peak$gc_content,
    nbins = nbins,
    show = TRUE,
    main = 'Fitted density for consensus peaks'
  )

  #  Do interpolation
  message('Interpolating for control regions ...')

  x_breaks = fit$x.breaks
  y_breaks = fit$y.breaks
  x_centers = (x_breaks[-1] + x_breaks[-length(x_breaks)]) / 2
  y_centers = (y_breaks[-1] + y_breaks[-length(y_breaks)]) / 2

  density_estimate = fields::interp.surface(
    obj = list(x = x_centers, y = y_centers, z = fit$counts),
    loc = cbind(log10(1 + meta_control_region$dist), meta_control_region$gc_content)
  )
  density_estimate[is.na(density_estimate)] = 0
  message(paste0(sum(density_estimate > 0), ' controls peaks with a nonzero weight.'))

  # Sampling
  message('Doing weighted random sampling ...')
  ipeaks = sample(nrow(meta_control_region), size = n_control, prob = density_estimate, replace = T)

  return(meta_control_region[ipeaks, ])
}
