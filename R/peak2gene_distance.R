#' Calculate Peak to Gene Distance
#'
#' This function calculates the distance between the center of peaks and the nearest gene transcription start site (TSS)
#' if the distance is less than a specified maximum distance.
#'
#' @param peak A tibble containing peak data with columns: 'chr', 'start', and 'end'.
#' @param ncores An integer specifying the number of cores to use for parallel processing.
#' @param meta_gene A tibble containing gene metadata with columns: 'id', 'symbol', 'biotype', 'chr', 'start', 'end', 'strand', 'tss' (default is grch38).
#' @param max_distance A numeric value specifying the maximum distance to consider (default is 500,000).
#' @param verbose A logical value indicating whether to print progress messages (default is TRUE).
#'
#' @return A tibble with columns: 'gene', 'peak', and 'distance'.
#' @export
#'
calc_peak2gene_distance <- function(peak, ncores, meta_gene = grch38, max_distance = 500e3, verbose = TRUE) {
  # Calculate peak center & name
  peak = peak %>%
    dplyr::mutate(center = (start + end) / 2) %>%
    dplyr::mutate(name = paste0('chr', chr, ':', start, '-', end))

  # Calculate the distance between centers of peaks & nearest gene TSS if < max_distance
  df = NULL
  for (chr in unique(meta_gene$chr)) {
    if (verbose) message(paste0('Working on chr ', chr, '...'))

    meta_gene_chr = meta_gene[meta_gene$chr == chr, ]
    peak_chr = peak[peak$chr == chr, ]

    df_chr = do.call(rbind, parallel::mcmapply(function(gene) {
      curr_meta_gene = meta_gene_chr %>% dplyr::filter(id == gene)
      res = NULL
      for (tss in curr_meta_gene$tss) {
        distances = abs(peak_chr$center - tss)
        selected = distances < max_distance
        if (any(selected)) {
          res = rbind(res, dplyr::tibble(peak = peak_chr$name[selected], distance = distances[selected]))
        }
      }
      if (is.null(res)) {
        return(NULL)
      } else {
        res = res %>%
          dplyr::group_by(peak) %>%
          dplyr::summarise(distance = min(distance), .groups = 'drop')
        return(dplyr::tibble(gene = gene, res))
      }
    }, unique(meta_gene_chr$id), SIMPLIFY = F, mc.cores = ncores))

    df = rbind(df, dplyr::tibble(chr = chr, df_chr))
  }

  return(df)
}
