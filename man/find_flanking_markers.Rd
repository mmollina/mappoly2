% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{find_flanking_markers}
\alias{find_flanking_markers}
\title{Find Flanking Markers}
\usage{
find_flanking_markers(A, B1, B2)
}
\arguments{
\item{A}{A character vector; positions of elements in this vector are considered.}

\item{B1}{A character vector; contains reference elements.}

\item{B2}{A character vector; elements for which closest flanking elements in \code{B1} are to be found.}
}
\value{
A list where each element in \code{B2} is mapped to its closest preceding and succeeding elements from \code{B1}.
}
\description{
Identifies the closest preceding and succeeding elements from \code{B1} for each element in \code{B2},
based on their positions in \code{A}.
}
\keyword{internal}
