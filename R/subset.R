#' Subset a mappoly2.data Object
#'
#' Creates a subset of individuals or markers from a \code{mappoly2.data} object based on specified criteria.
#'
#' @param x A \code{mappoly2.data} object from which the subset will be created.
#' @param type Character string specifying whether to subset "marker" or "individual". Defaults to "marker".
#' @param perc Numeric value indicating the proportion of individuals or markers to sample if \code{n} is \code{NULL}. Defaults to 0.1.
#' @param n Integer specifying the number of individuals or markers to sample. If \code{NULL}, \code{perc} must be specified.
#' @param select.ind Character vector of individual names to select. Used only if \code{type = "individual"}, and both \code{n} and \code{perc} are \code{NULL}.
#' @param select.mrk Character vector of marker names to select. Used only if \code{type = "marker"}, and both \code{n} and \code{perc} are \code{NULL}.
#' @param seed Integer or \code{NULL}; sets the seed for random number generation for reproducible sampling.
#' @param filter.non.conforming Logical; if \code{TRUE} (default), data points with unexpected genotypes (e.g., double reduction) are converted to \code{NA}.
#' @param filter.redundant Logical; if \code{TRUE} (default), removes redundant markers during map construction, keeping them annotated for export in the final map.
#' @param ... Additional arguments (currently not used).
#'
#' @return A \code{mappoly2.data} object containing the selected subset of individuals or markers.
#'
#' @details This function allows for flexible subsetting of \code{mappoly2.data} objects based on the number or proportion of individuals or markers. It supports random sampling and specific selection based on provided lists. Additional filtering options are available to handle non-conforming data and redundant markers.
#'
#' @examples
#' \donttest{
#' # Subset 20% of markers
#' subset_data <- subset.mappoly2.data(mappoly_data, type = "marker", perc = 0.2)
#'
#' # Subset 50 individuals
#' subset_data <- subset.mappoly2.data(mappoly_data, type = "individual", n = 50)
#'
#' # Select specific markers
#' subset_data <- subset.mappoly2.data(mappoly_data, select.mrk = c("marker1", "marker2"))
#' }
#'
#' @importFrom assertthat assert_that
#' @export
subset.mappoly2.data <- function(x,
                                 type = c("marker", "individual"),
                                 perc = 0.1,
                                 n = NULL,
                                 select.mrk = NULL,
                                 select.ind = NULL,
                                 seed = NULL,
                                 filter.non.conforming = TRUE,
                                 filter.redundant = TRUE,
                                 ...) {
  assert_that(is.mappoly2.data(x))

  # Determine sub-setting type based on input
  if (!is.null(select.ind)) {
    type <- "individual"
  } else {
    type <- match.arg(type)
  }

  # Initialize the subset object
  y <- x

  if (type == "marker") {
    # Subset markers
    if (is.null(select.mrk)) {
      if (!is.null(seed)) set.seed(seed)
      total_markers <- length(x$mrk.names)
      if (is.null(n)) {
        n <- ceiling(total_markers * perc)
      }
      if (n >= total_markers) {
        message("Number of selected markers is greater than or equal to the total number of markers. Returning original dataset.")
      } else {
        select.mrk <- sample(x$mrk.names, n)
        select.mrk <- select.mrk[na.omit(match(x$mrk.names, select.mrk))]
        y <- subset_data(x, select.mrk = select.mrk,
                         filter.non.conforming = filter.non.conforming,
                         filter.redundant = filter.redundant)
      }
    } else {
      y <- subset_data(x, select.mrk = select.mrk,
                       filter.non.conforming = filter.non.conforming,
                       filter.redundant = filter.redundant)
    }
  } else if (type == "individual") {
    # Subset individuals
    if (is.null(select.ind)) {
      if (!is.null(seed)) set.seed(seed)
      total_individuals <- length(x$ind.names)
      if (is.null(n)) {
        n <- ceiling(total_individuals * perc)
      }
      if (n >= total_individuals) {
        message("Number of selected individuals is greater than or equal to the total number of individuals. Returning original dataset.")
      } else {
        select.ind <- sample(x$ind.names, n)
        select.ind <- select.ind[na.omit(match(x$ind.names, select.ind))]
        y <- subset_data(x, select.ind = select.ind,
                         filter.non.conforming = filter.non.conforming,
                         filter.redundant = filter.redundant)
      }
    } else {
      y <- subset_data(x, select.ind = select.ind,
                       filter.non.conforming = filter.non.conforming,
                       filter.redundant = filter.redundant)
    }
  }

  # Create metadata for the subsetting step
  map_step <- data.frame(
    Metric = c(
      "Input markers",
      "Output markers",
      "Input individuals",
      "Output individuals",
      "Percentage of missing in output"
    ),
    Value = c(
      x$n.mrk,
      y$n.mrk,
      x$n.ind,
      y$n.ind,
      paste0(round(100 * sum(is.na(y$geno.dose)) / length(y$geno.dose), 1), "%")
    ),
    stringsAsFactors = FALSE
  )

  # Update the metadata of the mappoly object
  y <- update_metadata(y,
                       map_step = map_step,
                       class_suffix = NULL)

  return(y)
}

subset_data <- function(x,
                        select.ind = NULL,
                        select.mrk = NULL,
                        filter.non.conforming = TRUE,
                        filter.redundant = TRUE) {
  assert_that(!is.null(select.ind) || !is.null(select.mrk))

  # Validate and set selected individuals
  if (is.null(select.ind)) {
    select.ind <- x$ind.names
  } else {
    assert_that(all(select.ind %in% x$ind.names))
  }

  # Validate and set selected markers
  if (is.null(select.mrk)) {
    select.mrk <- x$mrk.names
  } else {
    assert_that(all(select.mrk %in% x$mrk.names))
  }

  # Create a copy of the original object to modify
  res <- x

  # Update the object with selected markers and individuals
  res$n.ind <- length(select.ind)
  res$n.mrk <- length(select.mrk)
  res$ind.names <- select.ind
  res$mrk.names <- select.mrk
  res$dosage.p1 <- res$dosage.p1[select.mrk]
  res$dosage.p2 <- res$dosage.p2[select.mrk]
  res$chrom <- res$chrom[select.mrk]
  res$genome.pos <- res$genome.pos[select.mrk]
  res$ref <- res$ref[select.mrk]
  res$alt <- res$alt[select.mrk]
  res$all.mrk.depth <- res$all.mrk.depth[select.mrk]
  res$geno.dose <- res$geno.dose[select.mrk, select.ind, drop = FALSE]

  # Screen non-conforming markers if required
  if (filter.non.conforming) {
    res <- filter_non_conforming_classes(res)
  }

  # Screen redundant markers if required
  if (filter.redundant) {
    res <- suppressMessages(filter_redundant(res, plot = FALSE))
  }

  # Update QAQC values
  read.depth <- NA
  if (!is.null(x$QAQC.values$markers)) {
    if ("read.depth" %in% colnames(x$QAQC.values$markers)) {
      read.depth <- x$QAQC.values$markers[res$mrk.names, "read.depth", drop = FALSE]
    }
  }

  full.sib <- NA
  if (!is.null(x$QAQC.values$individuals)) {
    if ("full.sib" %in% colnames(x$QAQC.values$individuals)) {
      full.sib <- x$QAQC.values$individuals[res$ind.names, "full.sib", drop = FALSE]
    }
  }

  res$QAQC.values <- .setQAQC(
    id.mrk = res$mrk.names,
    id.ind = res$ind.names,
    miss.mrk = rowMeans(is.na(res$geno.dose)),
    miss.ind = colMeans(is.na(res$geno.dose)),
    chisq.pval = suppressWarnings(mappoly_chisq_test(res)),
    read.depth = read.depth,
    full.sib = full.sib
  )

  return(res)
}
