% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_buster_score.R
\name{read_null_cbscore}
\alias{read_null_cbscore}
\title{Read Null Cluster-Buster Scores}
\usage{
read_null_cbscore(motif, dir_null_cbscore)
}
\arguments{
\item{motif}{A character string representing the motif name.}

\item{dir_null_cbscore}{A character string representing the directory containing the null cluster-buster scores.}
}
\value{
A named numeric vector of scores.
}
\description{
This function reads the cluster-buster scores on control regions for a given motif.
}
