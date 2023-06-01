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
  if (is.mappoly2.sequence(x)) {
    if(all(x$data$dosage.p2 == 0 | x$data$dosage.p2 == x$data$ploidy.p2))
      return("p1")
    if(all(x$data$dosage.p1 == 0 | x$data$dosage.p1 == x$data$ploidy.p1))
      return("p2")
    else
      return("both")
  }
  # else if (inherits(x, "mappoly.map")) {
  #   if(all(x$info$seq.dose.p2 == 0 | x$info$seq.dose.p2 == x$info$ploidy))
  #     return("p1")
  #   else if(all(x$info$seq.dose.p1 == 0 | x$info$seq.dose.p1 == x$info$ploidy))
  #     return("p2")
  #   else
  #     return("both")
  # }
  else{
    stop(deparse(substitute(x)), " is not an object of class 'mappoly2.sequence'")
  }
}
