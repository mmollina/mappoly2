get_seq_indices <- function(input.seq){
  match(input.seq$mrk.names, input.seq$data$mrk.names)
}

get_QAQCmrk_indices <- function(x){
  match(x$screened.data$mrk.names, x$mrk.names)
}

get_mrk_indices_from_chrom <- function(x, chrom){
  assert_that(mappoly2:::has.chromosome.info(x))
  ch.n.arg <- mappoly2:::embedded_to_numeric(chrom)
  assert_that(!any(is.na(ch.n.arg)), msg = "provide a valid chromosome identifier")
  ch.n.dat <- mappoly2:::embedded_to_numeric(x$chrom)
  ch.id <- intersect(which(ch.n.dat%in%ch.n.arg), mappoly2:::get_QAQCmrk_indices(x))
  ch.id
}


embedded_to_numeric <- function(x) {
  as.integer(gsub("[^0-9]", "", x))
}

#' Msg function
#'
#' @param void internal function to be documented
#' @keywords internal
#' @importFrom cli rule
#' @importFrom magrittr %>%
msg <- function(text, line = 1, col = "-"){
  cli::rule(line = line, left = text) %>%
    text_col(col = col) %>%
    message()
}

#' @importFrom rstudioapi isAvailable hasFun getThemeInfo
#' @importFrom crayon white black
text_col <- function(x, col = c("-","red", "blue")) {
  col <- match.arg(col)
  if(col == "red") return(crayon::red(x))
  else if(col == "blue") return(crayon::blue(x))
  # If RStudio not available, messages already printed in black
  if (!rstudioapi::isAvailable()) {
    return(x)
  }
  if (!rstudioapi::hasFun("getThemeInfo")) {
    return(x)
  }
  theme <- rstudioapi::getThemeInfo()
  if (isTRUE(theme$dark)) crayon::white(x) else crayon::black(x)
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

print_matrix <- function(mat, spaces = 5, zero.print = ".", row.names = TRUE){
  mat[mat==0]<-zero.print
  txt1 <- NULL
  for(i in 1:ncol(mat))
    txt1 <- c(txt1, colnames(mat)[i], mat[,i])
  n1 <- nchar(txt1)
  for (i in 1:length(txt1))
    txt1[i] <- paste(txt1[i], paste0(rep(" ", max(n1) - n1[i]), collapse = ""))
  dim(txt1) <- c(nrow(mat)+1, ncol(mat))
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
  for(i in 2:nrow(txt1)){
    if(row.names) cat(txt2[i], txt1[i,], "\n")
    else cat(paste0(rep(" ", nchar(txt2[i])), collapse = ""), txt1[i,], "\n")
  }
  cat(paste0(paste0(rep(" ", spaces), collapse = ""),
             paste0(rep("-", sum(nchar(c(txt2[1], txt1[1,])))+ncol(mat)-spaces), collapse = "")), "\n")
}

get_markers_from_chromosome <- function(x, arg){
  assert_that(has.chromosome.info(x))
  pattern <- "(ch|chr|CH|Chr|CHR|chrom|Chrom|Chromsome)"
  if(is.numeric(arg))
    ch.n.arg <- arg
  else{
    assert_that(all(is.character(arg)))
    assert_that(sum(grepl(pattern, arg, ignore.case = TRUE))  ==  length(arg))
    ch.n.arg <- embedded_to_numeric(arg)
  }
  ch.n.dat <- embedded_to_numeric(x$data$chrom)
  ch.id <- ch.n.dat%in%ch.n.arg
  mrk.names <- x$mrk.names[ch.id]
  chrom <- x$data$chrom[mrk.names]
  data.frame(chrom)
}

get_markers_from_grouped_sequence <- function(x, arg){
  assert_that(is.grouped.sequence(x))
  assert_that(is.numeric(arg))
  names(x$linkage.groups$groups.snp)[x$linkage.groups$groups.snp  %in%  arg]
}

get_markers_from_grouped_and_chromosome <- function(x, lg, ch = NULL){
  assert_that(is.pairwise.sequence(x))
  mrk.id.ch <- x$mrk.names
  if(!is.null(ch))
    mrk.id.ch <- rownames(mappoly2:::get_markers_from_chromosome(x, ch))
  assert_that(length(mrk.id.ch) > 0)
  mrk.id.lg <- mappoly2:::get_markers_from_grouped_sequence(x, lg)
  return(intersect(mrk.id.ch, mrk.id.lg))
}




