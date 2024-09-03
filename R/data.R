#' GRCh38 gene meta data
#'
#' Genes on chromosome 1-22, X, Y. Extracted from biomaRt on Sep 3, 2024. See data.R for details.
#'
"grch38"

# ensembl = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#
# grch38 = biomaRt::getBM(
#   attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype",
#                  "chromosome_name", "start_position", "end_position",
#                  "strand", "transcription_start_site"),
#   mart = ensembl
# )
#
# grch38 = dplyr::as_tibble(meta_gene) %>%
#   dplyr::rename(
#     id = ensembl_gene_id, symbol = external_gene_name, biotype = gene_biotype,
#     chr = chromosome_name, start = start_position, end = end_position,
#     tss = transcription_start_site
#   ) %>%
#   dplyr::filter(chr %in% c(as.character(1:22), 'X', 'Y'))
