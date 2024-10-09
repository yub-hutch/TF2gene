#' Read Narrow Peaks File
#'
#' This function reads a narrow peaks file into a tibble.
#'
#' @param file A string specifying the path to the narrow peaks file.
#' @return A tibble containing the narrow peaks data.
#' @export
read_narrow_peaks <- function(file) {
  x = data.table::fread(file, header = F)
  names(x) = c("Chromosome", "Start", "End", "Name", "Score", "Strand",
               "FC_summit", "-log10_pval", "-log10_qval", "Summit")
  x = dplyr::as_tibble(x)
  return(x)
}


#' Write Narrow Peaks File
#'
#' This function writes a tibble to a narrow peaks file.
#'
#' @param x A tibble containing the narrow peaks data.
#' @param file A string specifying the path to the output file.
#' @return None
#' @export
write_narrow_peaks <- function(x, file) {
  write.table(x, file = file, sep = '\t', quote = F, row.names = F, col.names = F)
}
