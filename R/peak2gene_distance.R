#' Calculate Peak to Gene Distance
#'
#' This function calculates the distance between the center of peaks and the nearest gene transcription start site (TSS)
#' if the distance is less than a specified maximum distance.
#'
#' @param peak A tibble containing peak data with columns: 'chr', 'start', 'end', and 'openness'.
#' @param meta_gene A tibble containing gene metadata with columns: 'id', 'symbol', 'biotype', 'chr', 'start', 'end', 'strand', 'tss'. For example, grch38.
#' @param max_distance A numeric value specifying the maximum distance to consider (default is 500,000).
#' @param ncores An integer specifying the number of cores to use for parallel processing.
#' @param verbose A logical value indicating whether to print progress messages.
#'
#' @return A tibble with columns: 'gene', 'peak', and 'distance'.
#' @export
#'
#' @examples
#' peak_data <- tibble::tibble(
#'   chr = c("1", "1", "X"),
#'   start = c(100, 200, 300),
#'   end = c(150, 250, 350),
#'   openness = c(0.5, 0.6, 0.7)
#' )
#'
#' gene_data <- tibble::tibble(
#'   id = c("gene1", "gene2", "gene3"),
#'   symbol = c("G1", "G2", "G3"),
#'   biotype = c("protein_coding", "protein_coding", "protein_coding"),
#'   chr = c("1", "1", "X"),
#'   start = c(50, 150, 250),
#'   end = c(200, 300, 400),
#'   strand = c("+", "-", "+"),
#'   tss = c(75, 175, 275)
#' )
#'
#' distances <- calc_peak2gene_distance(peak_data, gene_data, max_distance = 500e3, ncores = 2, verbose = TRUE)
calc_peak2gene_distance <- function(peak, meta_gene, max_distance = 500e3, ncores, verbose) {
  # Calculate peak center & name
  peak = peak %>%
    dplyr::mutate(center = (start + end) / 2) %>%
    dplyr::mutate(name = paste0('chr', chr, ':', start, '-', end))

  # Calculate the distance between peaks & nearest gene TSS if < max_distance
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

    df = rbind(df, df_chr)
  }

  return(df)
}
