#' Multipoint analysis using Hidden Markov Models
#'
#' @param void internal function
#' @keywords internal
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
hmm_map_reconstruction <- function(x,
                                   p,
                                   rf = NULL,
                                   verbose = FALSE,
                                   tol = 10e-4)
{
  assert_that(all(length(x$mrk.names) == sapply(p, nrow)))
  assert_that(x$data$ploidy.p1 == ncol(p[[1]]))
  assert_that(x$data$ploidy.p2 == ncol(p[[2]]))
  assert_that(all(apply(p[[1]], 1, sum) == x$data$dosage.p1[x$mrk.names]))
  assert_that(all(apply(p[[2]], 1, sum) == x$data$dosage.p2[x$mrk.names]))
  g <- x$data$geno.dose[x$mrk.names, ]
  pedigree <- matrix(rep(c(1,
                           2,
                           x$data$ploidy.p1,
                           x$data$ploidy.p2, 1),
                         x$data$n.ind),
                     nrow = x$data$n.ind,
                     byrow = TRUE)
  ret_H0 <- TRUE
  if(is.null(rf)){
    rf <- rep(0.01, nrow(g) - 1)
    ret_H0 <- FALSE
  }
  assert_that(length(rf) == nrow(g) - 1)
  w <- est_hmm_map_biallelic(PH = p,
                               G = g,
                               pedigree = pedigree,
                               rf = rf,
                               verbose = verbose,
                               tol = tol,
                               ret_H0 = ret_H0)
  names(w) <- c("loglike", "rec.frac")
  return(w)
}