% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bed9_importer.R
\name{bed9_importer}
\alias{bed9_importer}
\title{Imports bed9 (9 column) files into R}
\usage{
bed9_importer(bed, uncompress_keep = TRUE, skip = 0)
}
\arguments{
\item{bed}{Path to the bed file}

\item{uncompress_keep}{If the bed file is compressed, should the uncompressed file be
retained or removed? Defaults to TRUE.}

\item{skip}{How many lines to skip? Defaults to zero. You must look at the
first few lines of the bed file to determine if non bed entries are present.}
}
\description{
Imports bed9 (9 column) files into R
}
\details{
\code{bed9_importer} imports a 9 column bed file and retains proper
encoding (character, integer) of the columns.
}
