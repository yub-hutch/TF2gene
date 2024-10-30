#' Calculate the distance between consensus peaks to genes
#'
#' This function calculates the distance between the center of peaks and the nearest TSS,
#' if the distance is less than a specified maximum distance.
#'
#' @param consensus_peak A tibble containing consensus peak data with columns: 'chr', 'start', and 'end'.
#' @param ncores Number of cores to use for parallel processing.
#' @param meta_gene A tibble containing gene metadata with columns: 'id', 'symbol', 'biotype', 'chr', 'start', 'end', 'strand', 'tss' (default is grch38).
#' @param max_distance Maximum distance to consider (default is 500,000).
#'
#' @return A sparse matrix representing the distances. To distinguish distance 0 from empty sparse matrix
#' elements, distance 0 is set to 0.1.
#' @export
#'
calc_peak2gene_distance <- function(consensus_peak, ncores, meta_gene = grch38, max_distance = 500e3) {
  # Calculate peak center & name
  message('Formatting inputs ...')
  consensus_peak = consensus_peak %>%
    dplyr::mutate(center = (start + end) / 2) %>%
    dplyr::mutate(name = paste0(chr, ':', start, '-', end))

  # Calculate the distance between peak centers & nearest TSS if < max_distance
  df = NULL
  for (chr in unique(consensus_peak$chr)) {
    message(paste0('Working on ', chr, '...'))
    meta_gene_chr = meta_gene[meta_gene$chr == strsplit(chr, 'chr')[[1]][2], ]
    consensus_peak_chr = consensus_peak[consensus_peak$chr == chr, ]
    df_chr = do.call(rbind, parallel::mcmapply(function(gene) {
      curr_meta_gene = meta_gene_chr %>% dplyr::filter(id == gene)
      res = NULL
      for (tss in curr_meta_gene$tss) {
        distances = abs(consensus_peak_chr$center - tss)
        selected = distances < max_distance
        if (any(selected)) {
          res = rbind(res, dplyr::tibble(
            consensus_peak = consensus_peak_chr$name[selected],
            distance = distances[selected]
          ))
        }
      }
      if (is.null(res)) return(NULL)
      res %>%
        dplyr::group_by(consensus_peak) %>%
        dplyr::summarise(distance = min(distance), .groups = 'drop') %>%
        dplyr::mutate(gene = gene)
    }, unique(meta_gene_chr$id), SIMPLIFY = F, mc.cores = ncores))
    df = rbind(df, dplyr::tibble(chr = chr, df_chr))
  }

  # Create sparse matrix
  message('Creating sparse matrix ...')
  consensus_peaks = consensus_peak$name
  genes = sort(unique(meta_gene$id))
  mat = Matrix::sparseMatrix(
    i = match(df$consensus_peak, consensus_peaks),
    j = match(df$gene, genes),
    x = pmax(df$distance, 0.1),
    dims = c(length(consensus_peaks), length(genes)),
    dimnames = list(consensus_peaks, genes)
  )

  return(mat)
}
