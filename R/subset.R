#' Subsets mappoly2.data Object
#'
#' This function selects a subset of either individuals or markers from an object of the `mappoly2.data` class.
#'
#' @param x A `mappoly2.data` object from which the subset is to be selected.
#' @param n An integer specifying the number of individuals or markers to be sampled. If NULL, `percentage` must be specified.
#' @param percentage A numeric value between 0 and 1 indicating the proportion of individuals or markers to be sampled.
#' This parameter is considered only if `n` is set to NULL.
#' @param type A character string, either "marker" or "individual", indicating whether markers or individuals should be subset.
#' @param selected.ind A character vector containing the names of the individuals to select.
#' This parameter takes effect only if `type = "individual"`, and both `n` and `percentage` are NULL.
#' @param selected.mrk A character vector containing the names of the markers to select.
#' This parameter takes effect only if `type = "marker"`, and both `n` and `percentage` are NULL.
#'
#' @return A `mappoly2.data` object that contains the selected subset of individuals or markers.
#' @export
subset.mappoly2.data <- function(x, type = c("marker", "individual"),
                                 perc = 0.1, n = NULL,
                                 select.mrk = NULL, select.ind = NULL,
                                 seed = NULL){
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
    return(subset_data(x, select.mrk = select.mrk))
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
    return(subset_data(x, select.ind = select.ind))
  }
}

subset_data <- function(input.data,
                        select.ind = NULL,
                        select.mrk = NULL){
  assert_that(!is.null(select.ind) | !is.null(select.mrk))
  assert_that(all(select.mrk%in%input.data$mrk.names))
  assert_that(all(select.ind%in%input.data$ind.names))
  if(is.null(select.ind))
    select.ind <- input.data$ind.names
  if(is.null(select.mrk))
    select.mrk <- input.data$mrk.names
  input.data$n.ind <- length(select.ind)
  input.data$n.mrk <- length(select.mrk)
  input.data$ind.names <- select.ind
  input.data$mrk.names <- select.mrk
  input.data$dosage.p1 <- input.data$dosage.p1[select.mrk]
  input.data$dosage.p2 <- input.data$dosage.p2[select.mrk]
  input.data$chrom <- input.data$chrom[select.mrk]
  input.data$genome.pos <- input.data$genome.pos[select.mrk]
  input.data$ref <- input.data$ref[select.mrk]
  input.data$alt <- input.data$alt[select.mrk]
  input.data$all.mrk.depth <- input.data$all.mrk.depth[select.mrk]
  input.data$geno.dose <- input.data$geno.dose[select.mrk,select.ind]
  return(input.data)
}
