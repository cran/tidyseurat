% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{get_abundance_sc_wide}
\alias{get_abundance_sc_wide}
\title{get abundance long}
\usage{
get_abundance_sc_wide(
  .data,
  transcripts = NULL,
  all = FALSE,
  assay = .data@active.assay,
  slot = "data"
)
}
\arguments{
\item{.data}{A tidyseurat}

\item{transcripts}{A character}

\item{all}{A boolean}

\item{...}{Parameters to pass to join wide, i.e. assay name to extract transcript abundance from}
}
\value{
A Seurat object
}
\description{
get abundance long
}
\keyword{internal}
