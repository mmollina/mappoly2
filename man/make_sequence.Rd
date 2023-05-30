% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sequence.R
\name{make_sequence}
\alias{make_sequence}
\alias{print.mappoly2.sequence}
\alias{plot.mappoly2.sequence}
\title{Create a sequence of markers}
\usage{
make_sequence(
  input.obj,
  arg = NULL,
  info.parent = c("all", "p1", "p2"),
  genomic.info = NULL,
  phase = NULL
)

\method{print}{mappoly2.sequence}(x, ...)

\method{plot}{mappoly2.sequence}(x, thresh.line = NULL, ...)
}
\arguments{
\item{input.obj}{an object of one of the following classes:
\code{mappoly.data}, \code{mappoly.map}, \code{mappoly.sequence},
\code{mappoly.group},
\code{mappoly.pcmap}, \code{mappoly.pcmap3d}, or \code{mappoly.geno.ord}}

\item{arg}{one of the following input types:
\enumerate{
\item A string 'all': Generates a sequence with all markers from the raw data.
\item A string or a vector of strings 'chrX': Specifies a chromosome (with 'X'
being the chromosome number). For unassigned markers, use 'chr0'.
\item A vector of integers: Indicates the position of markers in the original
data set to be included in the sequence.
\item A vector of strings: Indicates the names or identifiers of the genetic
markers
\item An integer: Represents a linkage group when input.object is of class
\code{mappoly.group}.
\item NULL: Applicable when \code{input.object} belongs to one of the following
classes: \code{mappoly.pcmap}, \code{mappoly.pcmap3d},
\code{mappoly.unique.seq}, or \code{mappoly.geno.ord}.
}}

\item{info.parent}{one of the following options:
\enumerate{
\item \code{'all'}{select all dosage combinations in both parents (default)}
\item \code{'P1'}{select informative markers parent 1}
\item \code{'P2'}{select informative markers parent 2}
}}

\item{genomic.info}{An optional argument applicable only to \code{mappoly.group}
objects, which can be either NULL or a numeric combination of sequences
from genomic information used to create the sequences:
\enumerate{
\item NULL (default): Returns a sequence containing all markers as defined by
the grouping function.
\item 1: Returns a sequence with markers that match the intersection between
the grouping function and genomic information, considering the genomic information
sequence with the maximum number of matching markers for the group.
\item c(1, 2): Returns a sequence with markers that match the intersection between
the grouping function and genomic information, considering the two genomic
information sequences with the maximum number of matching markers for the group,
and so on.
}}

\item{x}{an object of the class \code{mappoly.sequence}}

\item{...}{currently ignored}

\item{thresh.line}{position of a threshold line for p values of the segregation test (default = \code{0.05/n.mrk})}
}
\value{
An object of class \code{mappoly.sequence}, which is a
list containing the following components:
\item{seq.phases}{a \code{list} with the linkage phases between
markers in the sequence, in corresponding positions. \code{-1}
means that there are no defined linkage phases}
\item{seq.rf}{a \code{vector} with the recombination
frequencies between markers in the sequence. \code{-1} means
that there are no estimated recombination frequencies}
\item{loglike}{log-likelihood of the corresponding linkage
map}
\item{data.name}{name of the object of class
\code{mappoly.data} with the raw data}
\item{twopt}{name of the object of class \code{mappoly.twopt}
with the 2-point analyses. \code{-1} means that the twopt
estimates were not computed}
}
\description{
Makes a sequence of markers based on an object of another class.
}
\author{
Marcelo Mollinari (\email{mmollin@ncsu.edu}) and
Gabriel Gesteira, (\email{gdesiqu@ncsu.edu})
}