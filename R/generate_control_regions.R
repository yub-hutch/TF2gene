#' Generate Genome-Wide Control Regions
#'
#' This function generates control regions encompassing a given set of genes. It creates intervals around gene transcription start sites (TSS), resizes them, and splits them into non-overlapping windows of a specified length.
#'
#' @param meta_gene A tibble containing gene metadata with columns: 'id', 'symbol', 'biotype', 'chr', 'start', 'end', 'strand', 'tss' (default is grch38).
#' @param chr_lens A tibble containing chromosome lengths with columns: 'chr' and 'length' (default is grch38_chr_lens).
#' @param region_length A numeric value specifying the length of each control_region (default is 500).
#' @param radius A numeric value specifying the radius around the TSS to consider for generating control regions (default is 500,000).
#'
#' @return A tibble with columns: 'chr', 'start', and 'end' representing the generated control regions.
#' @export
#'
generate_control_regions <- function(meta_gene = grch38, chr_lens = grch38_chr_lens, region_length = 500, radius = 500e3) {
  control_region = NULL
  for (chr in unique(meta_gene$chr)) {
    meta_gene_chr = meta_gene[meta_gene$chr == chr, ]
    chr_len = chr_lens$length[chr_lens$chr == chr]

    # Get intervals encompassing all genes
    ir = meta_gene_chr %>%
      dplyr::group_by(id) %>%
      dplyr::summarise(
        start = max(1, min(tss) - radius),
        end = min(chr_len, max(tss) + radius),
        .groups = 'drop'
      ) %>%
      dplyr::select(start, end)

    ir = IRanges::IRanges(start = ir$start, end = ir$end)

    # Merge overlapping intervals
    ir = IRanges::reduce(ir)

    # Resize the intervals
    new_widths = floor(ir@width / region_length) * region_length
    ir = IRanges::resize(ir, new_widths, fix = 'center')
    ir = ir[IRanges::start(ir) > 0 & IRanges::end(ir) <= chr_len]

    # Split intervals into control regions
    region_starts = unlist(lapply(seq_along(ir), function(i) {
      seq(IRanges::start(ir[i]), IRanges::end(ir[i]), by = region_length)
    }))
    region_ends = region_starts + region_length - 1

    control_region = rbind(control_region, dplyr::tibble(chr = chr, start = region_starts, end = region_ends))
  }

  control_region$start = control_region$start - 1 # Align to bed format

  return(control_region)
}


#' Write control regions to BED File
#'
#' This function converts a tibble of control regions into BED format and writes it to a specified file.
#'
#' @param control_region A tibble containing control region data with columns: 'chr', 'start', and 'end'.
#' @param fbed A character string specifying the file path to write the BED file.
#'
#' @return None. The function writes the BED file to the specified path.
#' @export
#'
write_control_regions_to_bed <- function(control_region, fbed) {
  control_region$chr = paste0('chr', control_region$chr)
  control_region$name = paste0('control', seq(nrow(control_region)))
  control_region$openness = '.'
  write.table(control_region, file = fbed, row.names = F, col.names = F, quote = F, sep = '\t')
}
