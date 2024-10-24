#' Summarize Features of Genomic Regions
#'
#' This function reads a FASTA file containing sequences of genomic regions and a metadata table of gene information.
#' It calculates the distance between the centers of the regions and the nearest TSS, and the A/C/G/T frequencies of each region.
#'
#' @param fasta A character string specifying the path to the FASTA file containing sequences of consensus peaks.
#' @param ncores An integer specifying the number of cores to use for parallel processing.
#' @param meta_gene A data frame containing metadata of genes, including chromosome and TSS information (default is grch38).
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
summarize_features_of_regions <- function(fasta, ncores, meta_gene = grch38) {
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
