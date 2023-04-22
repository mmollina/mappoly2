#' Subsetting mappoly2 data set
#'
#' @param x an object  of class \code{mappoly.data}
#' @param selected.ind a vector containing the name of the individuals to select.
#' @param selected.mrk a vector containing the name of the markers to select.
#'
#' @return an object  of class \code{mappoly.data}
#' @keywords internal
#' @export
subset_data <- function(input.data,
                        select.ind = NULL,
                        select.mrk = NULL){
  assert_that(!is.null(select.ind) | !is.null(select.mrk))
  assert_that(all(select.mrk%in%input.data$mrk.names))
  assert_that(all(select.ind%in%input.data$ind.names))
  flg <- is.null(select.ind)
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
  input.data$chisq.pval <- input.data$chisq.pval[select.mrk]
  if(flg)
    return(structure(suppressWarnings(mappoly_chisq_test(input.data)), class = class(input.data)))
  return(input.data)
}
