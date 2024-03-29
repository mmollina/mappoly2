% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_methods.R
\name{plot_map}
\alias{plot_map}
\title{Plot Genetic Map}
\usage{
plot_map(
  x,
  lg = 1,
  type = c("mds", "genome"),
  parent = c("p1p2", "p1", "p2"),
  left.lim = 0,
  right.lim = Inf,
  phase = TRUE,
  mrk.names = FALSE,
  plot.dose = TRUE,
  homolog.names.adj = 3,
  cex = 1,
  xlim = NULL,
  main = "",
  ...
)
}
\arguments{
\item{x}{An object representing genetic mapping data.}

\item{lg}{The linkage group to be visualized, default is the first linkage group (lg = 1).}

\item{type}{The type of genetic map to process, either "mds" or "genome".}

\item{parent}{Specifies which parent's data to use in the visualization.
Options are "p1p2" (both parents), "p1" (first parent), or "p2" (second parent).}

\item{left.lim}{The left limit for the plotting area, default is 0.}

\item{right.lim}{The right limit for the plotting area, default is Inf.}

\item{phase}{Logical; if TRUE, phases are included in the plot.}

\item{mrk.names}{Logical; if TRUE, marker names are displayed on the plot.}

\item{plot.dose}{Logical; if TRUE, doses are plotted.}

\item{homolog.names.adj}{Adjustment for homolog names in the plot.}

\item{cex}{Character expansion size for text in the plot.}

\item{xlim}{The limits for the x-axis. Can be set to NULL for automatic adjustment.}

\item{main}{The main title for the plot.}

\item{...}{Additional graphical parameters.}
}
\value{
The function does not return a value but generates a plot of the genetic map.
}
\description{
This function visualizes a genetic map for a specified linkage group. It supports various types of genetic maps and offers customization options for the display.
}
\details{
The function creates a detailed plot of a genetic map for a given linkage group.
It can display various features such as phases, marker names, and doses, and allows
for customization of the plot's appearance.
}
