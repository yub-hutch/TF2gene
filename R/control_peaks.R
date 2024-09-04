#' Generate Null Peaks
#'
#' This function generates control peaks encompassing a given set of genes. It creates intervals around gene transcription start sites (TSS), resizes them, and splits them into non-overlapping windows of a specified length.
#'
#' @param meta_gene A tibble containing gene metadata with columns: 'id', 'symbol', 'biotype', 'chr', 'start', 'end', 'strand', 'tss' (default is grch38).
#' @param chr_lens A tibble containing chromosome lengths with columns: 'chr' and 'length' (default is hg38_chr_lens).
#' @param peak_length A numeric value specifying the length of each peak (default is 500).
#' @param radius A numeric value specifying the radius around the TSS to consider for generating peaks (default is 500,000).
#'
#' @return A tibble with columns: 'chr', 'start', and 'end' representing the generated null peaks.
#' @export
#'
get_control_peaks <- function(meta_gene = grch38, chr_lens = hg38_chr_lens, peak_length = 500, radius = 500e3) {
  null_peak = NULL
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
    new_widths = floor(ir@width / peak_length) * peak_length
    ir = IRanges::resize(ir, new_widths, fix = 'center')
    ir = ir[IRanges::start(ir) > 0 & IRanges::end(ir) <= chr_len]

    # Split intervals into null peaks
    peak_starts = unlist(lapply(seq_along(ir), function(i) {
      seq(IRanges::start(ir[i]), IRanges::end(ir[i]), by = peak_length)
    }))
    peak_ends = peak_starts + peak_length - 1

    null_peak = rbind(null_peak, dplyr::tibble(chr = chr, start = peak_starts, end = peak_ends))
  }

  return(null_peak)
}


#' Write Control Peaks to BED File
#'
#' This function converts a tibble of control peaks into BED format and writes it to a specified file.
#'
#' @param peak A tibble containing peak data with columns: 'chr', 'start', and 'end'.
#' @param fbed A character string specifying the file path to write the BED file.
#'
#' @return None. The function writes the BED file to the specified path.
#' @export
#'
write_control_peak_to_bed <- function(peak, fbed) {
  peak$start = peak$start - 1 # bed file is 0-based
  peak$chr = paste0('chr', peak$chr)
  peak$name = paste0('null', seq(nrow(peak)))
  peak$openness = '.'
  write.table(peak, file = fbed, row.names = F, col.names = F, quote = F, sep = '\t')
}
