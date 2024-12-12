#' GRCh38 gene meta data
#'
#' Genes on chromosome 1-22, X, Y. Extracted from biomaRt on Sep 3, 2024. Documented in data-raw/.
"grch38"


#' Chromosome length of GRCh38 build
#'
#' Only retain chromosome 1-22, X, Y. Documented in data-raw/.
"grch38_chr_lens"


#' Metadata of motifs
#'
#' A tibble of the metadata of motifs.
"meta_motif"


#' TF-motif annotations from cisTarget
#'
#' A tibble annotating connections between TF and motifs.
"df_tf2motif"


#' Processed TF-motif annotations
#'
#' A sparse TF by motif matrix. 3: Direct; 2: By orthology; 1: By similarity.
"mat_tf2motif"


#' Meta data of genome-wide control regions
#'
#' A tibble with columns chr, start, end, dist (distance to nearest TSS), freq_A, freq_C, freq_G, freq_T.
"meta_control_region"
