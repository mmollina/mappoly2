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

print_matrix <- function(mat, spaces = 5, zero.print = "."){
  mat[mat==0]<-zero.print
  txt1 <- c(colnames(mat), mat)
  n1 <- c(nchar(colnames(mat)), nchar(mat))
  for (i in 1:length(txt1))
    txt1[i] <- paste(txt1[i], paste0(rep(" ", max(n1) - n1[i]), collapse = ""))
  dim(txt1) <- c(nrow(mat), ncol(mat)+1)
  txt1<-t(txt1)
  txt2 <- rownames(mat)
  n2 <- nchar(txt2)
  for (i in 1:length(txt2))
    txt2[i] <- paste(txt2[i], paste0(rep(" ", max(n2) - n2[i]), collapse = ""))
  txt2 <- c(paste0(rep(" ", (nchar(max(n2))+1)), collapse = ""), txt2)
  for (i in 1:length(txt2))
    txt2[i] <- paste0(paste0(rep(" ", spaces), collapse = ""), txt2[i])
  cat(txt2[1], txt1[1,], "\n")
  cat(paste0(paste0(rep(" ", spaces), collapse = ""),
             paste0(rep("-", sum(nchar(c(txt2[1], txt1[1,])))+ncol(mat)-spaces), collapse = "")), "\n")
  for(i in 2:nrow(txt1))
    cat(txt2[i], txt1[i,], "\n")
  cat(paste0(paste0(rep(" ", spaces), collapse = ""),
             paste0(rep("-", sum(nchar(c(txt2[1], txt1[1,])))+ncol(mat)-spaces), collapse = "")), "\n")
}

