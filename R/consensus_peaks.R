#' Format ATAC-seq consensus peak
#'
#' This function reads a BED file and formats the consensus peaks. It only retains peaks in chromosome 1-22, X, & Y.
#'
#' @param fbed A character string specifying the path to the BED file.
#'
#' @return A tibble with four columns:
#' \describe{
#'   \item{chr}{Chromosome identifier.}
#'   \item{start}{Start position of the peak (0-based, included).}
#'   \item{end}{End position of the peak (0-based, excluded).}
#'   \item{openness}{Openness score of the peak.}
#' }
#' @export
#'
format_consensus_peaks <- function(fbed) {
  peak = dplyr::as_tibble(data.table::fread(fbed, header = F))[, c(1, 2, 3, 5)]
  names(peak) = c('chr', 'start', 'end', 'openness')

  peak = peak %>%
    dplyr::filter(chr %in% paste0('chr', c(as.character(1:22), 'X', 'Y'))) %>%
    dplyr::mutate(chr = sapply(strsplit(chr, 'chr'), tail, n = 1),
                  name = paste0('chr', chr, ':', start, '-', end))

  return(peak)
}


#' Summarize ATAC-seq Consensus Peaks
#'
#' This function reads a FASTA file containing sequences of consensus peaks and a metadata table of gene information. It calculates the distance between the centers of the peaks and the nearest gene transcription start site (TSS), and computes the GC content of each peak.
#'
#' @param fasta A character string specifying the path to the FASTA file containing sequences of consensus peaks.
#' @param ncores An integer specifying the number of cores to use for parallel processing.
#' @param meta_gene A data frame containing metadata of genes, including chromosome and TSS information (default is grch38).
#' @param verbose A logical value indicating whether to print progress messages (default is TRUE).
#'
#' @return A tibble with the following columns:
#' \describe{
#'   \item{chr}{Chromosome identifier.}
#'   \item{start}{Start position of the peak (0-based, included).}
#'   \item{end}{End position of the peak (0-based, excluded).}
#'   \item{name}{Name of the peak.}
#'   \item{center}{Center position of the peak.}
#'   \item{dist}{Distance from the center of the peak to the nearest gene TSS.}
#'   \item{gc_content}{GC content of the peak sequence.}
#' }
#' @export
#'
summarize_consensus_peaks <- function(fasta, ncores, meta_gene = grch38, verbose = TRUE) {
  # Read in sequences of consensus peaks
  if (verbose) message('Loading peak sequences ...')
  fa = seqinr::read.fasta(fasta, seqtype = 'DNA', as.string = F, forceDNAtolower = T)

  # Get genomic locations of consensus peaks
  peak_names = names(fa)
  peak = dplyr::as_tibble(t(sapply(peak_names, function(peak_name) {
    tmp = strsplit(peak_name, '[:-]')[[1]]
    chr = strsplit(tmp[1], 'chr')[[1]][2]
    c(chr = chr, start = tmp[2], end = tmp[3])
  }))) %>% dplyr::mutate(
    name = peak_names,
    start = as.integer(start),
    end = as.integer(end),
    center = (start + end) / 2
  )

  # Calculate the distance between centers peaks & nearest gene TSS
  dist = NULL
  for (chr in unique(meta_gene$chr)) {
    if (verbose) message(paste0('Working on chr ', chr, ' ...'))
    meta_gene_chr = meta_gene[meta_gene$chr == chr, ]
    peak_chr = peak[peak$chr == chr, ]

    peak_chr$dist = parallel::mcmapply(function(center) {
      min(abs(center - meta_gene_chr$tss))
    }, peak_chr$center, mc.cores = ncores)

    dist = rbind(dist, peak_chr)
  }

  # GC content
  gc_contents = parallel::mcmapply(function(seq) mean(seq %in% c('c', 'g')), fa, mc.cores = ncores)
  dist$gc_content = gc_contents[dist$name]

  return(dist)
}


#' Plot Metadata of Peaks
#'
#' This function creates a hexbin plot of the distance to the nearest transcription start site (TSS) versus the GC content of the peak sequences.
#'
#' @param meta_peak A tibble containing metadata of peaks, including distance to the nearest TSS and GC content. See \code{\link{summarize_consensus_peaks}}.
#' @param log_dist A logical value indicating whether to log-transform the distance to the nearest TSS (default is TRUE).
#' @param add_cor A logical value indicating whether to add correlation test results to the plot (default is TRUE).
#'
#' @return A ggplot object representing the hexbin plot.
#' @export
#'
plot_meta_peak <- function(meta_peak, log_dist = TRUE, add_cor = TRUE) {
  if (log_dist) meta_peak$dist = meta_peak$dist + 1
  p = ggplot2::ggplot(meta_peak, ggplot2::aes(dist, gc_content)) +
    ggplot2::geom_hex(bins = 100) +
    ggplot2::labs(x = 'Distance to nearest TSS', y = 'GC content') +
    ggplot2::guides(fill = 'none') +
    ggpubr::theme_pubr()
  if (log_dist) p = p + ggplot2::scale_x_log10()
  if (add_cor) p = p + ggpubr::stat_cor(method = 'spearman', label.x.npc = 'center')
  return(p)
}
