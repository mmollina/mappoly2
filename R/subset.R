#' Subset mappoly2.data Object
#'
#' This function creates a subset of either individuals or markers from an object
#' of the `mappoly2.data` class, based on specified criteria.
#'
#' @param x A `mappoly2.data` object from which the subset is to be selected.
#' @param type Character string, either "marker" or "individual", indicating whether
#'             markers or individuals should be subset. Default is "marker".
#' @param perc Numeric; the proportion of individuals or markers to be sampled.
#'             This is used only if `n` is NULL. Default is 0.1.
#' @param n Integer; the number of individuals or markers to be sampled.
#'           If NULL, `perc` must be specified.
#' @param select.ind Character vector containing the names of the individuals to select.
#'                   Effective only if `type = "individual"`, and both `n` and `perc` are NULL.
#' @param select.mrk Character vector containing the names of the markers to select.
#'                   Effective only if `type = "marker"`, and both `n` and `perc` are NULL.
#' @param seed Integer or NULL; sets the seed for random number generation for reproducible sampling.
#' @param filter.non.conforming Logical; if TRUE (default), data points with unexpected genotypes
#'                              (e.g., double reduction) are converted to 'NA'.
#' @param filter.redundant Logical; if TRUE (default), removes redundant markers during map construction,
#'                         keeping them annotated to export to the final map.
#' @param ... Additional arguments, currently not used.
#'
#' @return A `mappoly2.data` object that contains the selected subset of individuals or markers.
#'
#' @details The function allows for flexible subsetting of `mappoly2.data` objects based on
#'          the number or proportion of individuals or markers. It supports random sampling
#'          and specific selection based on provided lists. Additional filtering options are
#'          available to handle non-conforming data and redundant markers.
#'
#'
#' @importFrom assertthat assert_that
#' @export
subset.mappoly2.data <- function(x, type = c("marker", "individual"),
                                 perc = 0.1, n = NULL,
                                 select.mrk = NULL, select.ind = NULL,
                                 seed = NULL, filter.non.conforming = TRUE,
                                 filter.redundant = TRUE, ...){
  assert_that(is.mappoly2.data(x))
  if(!is.null(select.ind))
    type <- "individual"
  else type <- match.arg(type)
  if(type  ==  "marker"){
    if(is.null(select.mrk)){
      if (!is.null(seed)) set.seed(seed)
      if(is.null(n)) n <- length(x$mrk.names) * perc
      if(n > length(x$mrk.names)){
        message("number of selected markers is larger\nthan number of markers in the original\ndataset. Returning original dataset")
        return(x)
      }
      select.mrk = sample(x$mrk.names, ceiling(n))
    }
    return(subset_data(x, select.mrk = select.mrk,
                       filter.non.conforming = filter.non.conforming,
                       filter.redundant = filter.redundant))
  } else if(type == "individual"){
    if(is.null(select.ind)){
      if (!is.null(seed)) set.seed(seed)
      if(is.null(n)) n <- length(x$ind.names) * perc
      if(n > length(x$ind.names)){
        message("number of selected individuals is larger\nthan number of individuals in the original\ndataset. Returning original dataset")
        return(x)
      }
      select.ind = sample(x$ind.names, ceiling(n))
    }
    return(subset_data(x, select.ind = select.ind,
                       filter.non.conforming = filter.non.conforming,
                       filter.redundant = filter.redundant))
  }
}

subset_data <- function(x,
                        select.ind = NULL,
                        select.mrk = NULL,
                        filter.non.conforming = TRUE,
                        filter.redundant = TRUE){
  assert_that(!is.null(select.ind) | !is.null(select.mrk))
  assert_that(all(select.mrk%in%x$mrk.names))
  assert_that(all(select.ind%in%x$ind.names))
  if(is.null(select.ind))
    select.ind <- x$ind.names
  if(is.null(select.mrk))
    select.mrk <- x$mrk.names
  x$n.ind <- length(select.ind)
  x$n.mrk <- length(select.mrk)
  x$ind.names <- select.ind
  x$mrk.names <- select.mrk
  x$dosage.p1 <- x$dosage.p1[select.mrk]
  x$dosage.p2 <- x$dosage.p2[select.mrk]
  x$chrom <- x$chrom[select.mrk]
  x$genome.pos <- x$genome.pos[select.mrk]
  x$ref <- x$ref[select.mrk]
  x$alt <- x$alt[select.mrk]
  x$all.mrk.depth <- x$all.mrk.depth[select.mrk]
  x$geno.dose <- x$geno.dose[select.mrk,select.ind]
  res <- x
  # Screening non-conforming markers
  if (filter.non.conforming) {
    res <- filter_non_conforming_classes(res)
  }
  if(filter.redundant){
    res <- filter_redundant(res)
  }
  res$QAQC.values <- .setQAQC(id.mrk = res$mrk.names,
                              id.ind = res$ind.names,
                              miss.mrk = apply(res$geno.dose, 1, function(x) sum(is.na(x)))/res$n.ind,
                              miss.ind = apply(res$geno.dose, 2, function(x) sum(is.na(x)))/res$n.mrk,
                              chisq.pval = suppressWarnings(mappoly_chisq_test(res)))
  return(res)
}
