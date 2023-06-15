get_seq_indices <- function(input.seq){
  match(input.seq$mrk.names, input.seq$data$mrk.names)
}
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

sort_phase <- function(x, only.best = TRUE){
  assert_that(is.mappoly2.sequence(x))
  ph <- x$phases[sapply(x$phases, function(x) !is.null(x$loglike))]
  if(only.best){
    x$phases <- ph[which.max(sapply(ph, function(x) x$loglike))]
    return(x)
  } else {
    x$phases <- ph[order(sapply(ph, function(x) x$loglike), decreasing = TRUE)]
    return(x)
  }
}

find_flanking_markers <- function(A, B1, B2) {
  A <- as.character(A)
  B1 <- as.character(B1)
  B2 <- as.character(B2)

  # Match B1 elements to their positions in A
  B1_positions <- match(B1, A)
  names(B1_positions) <- B1

  flanking_letters <- lapply(B2, function(b) {
    b_position <- match(b, A)
    preceding <- B1_positions[B1_positions < b_position]
    succeeding <- B1_positions[B1_positions > b_position]

    # Get the closest preceding and succeeding elements from B1
    preceding <- if (length(preceding) > 0) names(preceding)[which.max(preceding)] else NA
    succeeding <- if (length(succeeding) > 0) names(succeeding)[which.min(succeeding)] else NA

    list(preceding = preceding, succeeding = succeeding)
  })

  names(flanking_letters) <- B2
  flanking_letters
}

