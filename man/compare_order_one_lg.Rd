% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map.R
\name{compare_order_one_lg}
\alias{compare_order_one_lg}
\title{Compare Marker Orders for a Single Linkage Group}
\usage{
compare_order_one_lg(
  g.mds,
  g.gen,
  ph.mds,
  ph.gen,
  rf.mds,
  rf.gen,
  parent,
  ploidy.p1,
  ploidy.p2,
  dosage.p1,
  dosage.p2,
  error = 0,
  verbose = TRUE,
  tol = 0.01
)
}
\arguments{
\item{g.mds}{Genotype data for MDS.}

\item{g.gen}{Genotype data for genome method.}

\item{ph.mds}{Phase information for MDS.}

\item{ph.gen}{Phase information for genome method.}

\item{rf.mds}{Recombination fractions for MDS.}

\item{rf.gen}{Recombination fractions for genome method.}

\item{parent}{Parent phase to consider.}

\item{ploidy.p1}{Ploidy of parent 1.}

\item{ploidy.p2}{Ploidy of parent 2.}

\item{dosage.p1}{Dosage information for parent 1.}

\item{dosage.p2}{Dosage information for parent 2.å}

\item{error}{error parameter.}

\item{verbose}{If \code{TRUE}, prints messages.}

\item{tol}{Tolerance level.}
}
\value{
A matrix with comparison results, including number of markers, map lengths,
maximum gaps, average distances, and log likelihoods for both methods.
}
\description{
This internal function compares the marker orders for a single linkage group
using two different map construction methods: multidimensional scaling (MDS)
and genome order. It calculates the mapping for each method and compares
their log likelihoods to determine the best map order.
}
\keyword{internal}
