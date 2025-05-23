% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matched_control_regions.R
\name{extract_matched_control_regions}
\alias{extract_matched_control_regions}
\title{Extract Matched control regions}
\usage{
extract_matched_control_regions(
  meta_consensus_peak,
  meta_control_region,
  n,
  same_chr,
  ncores
)
}
\arguments{
\item{meta_consensus_peak}{A tibble containing metadata of consensus peaks. See \code{\link{summarize_features_of_regions}}.}

\item{meta_control_region}{A tibble containing metadata for control regions.}

\item{n}{Number of matched control regions to sample for each consensus peak.}

\item{same_chr}{Logical indicating whether to select control regions from the same chromosome (TRUE) or any chromosome (FALSE).}

\item{ncores}{An integer specifying the number of cores to use for parallel processing.}
}
\value{
A sparse matrix recording the distances between each consensus peak and its matched control regions.
Distance of 0 is replaced by 1e-10 for sparse matrix representation.
}
\description{
This function samples control regions that match the A/C/G/T frequencies of consensus peaks.
}
