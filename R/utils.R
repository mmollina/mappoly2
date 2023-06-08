get_seq_indices <- function(input.seq)
  match(input.seq$mrk.names, input.seq$data$mrk.names)

embedded_to_numeric <- function(x) {
  as.integer(gsub("[^0-9]", "", x))
}

#' Detects which parent is informative
#'
#' @param x an object of class \code{mappoly2.sequence}
#'
#' @export
detect_info_par<-function(x){
  ## checking for correct object
  assert_that(is.mappoly2.sequence(x))
    if(all(x$data$dosage.p2 == 0 | x$data$dosage.p2 == x$data$ploidy.p2))
      return("p1")
    if(all(x$data$dosage.p1 == 0 | x$data$dosage.p1 == x$data$ploidy.p1))
      return("p2")
    else
      return("both")
}

sort_phase <- function(x){
  assert_that(is.mappoly2.sequence(x))
  ph <- x$phases[sapply(x$phases, function(x) !is.null(x$loglike))]
  x$phases <- ph[order(sapply(ph, function(x) x$loglike), decreasing = TRUE)]
  return(x)
}
