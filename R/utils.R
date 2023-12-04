#' Given a screened data set, returns the indices of the
#' screened markers in the corresponding raw dataset
get_screened_mrk_indices <- function(x){
  match(x$screened.data$mrk.names, x$mrk.names)
}

#' Given a dataset containing chromosome information,
#' and a chromosome vector return a vector of screened
#' markers corresponding to the chromosomes provided
get_mrk_indices_from_chrom <- function(x, chrom){
  assert_that(mappoly2:::has.chromosome.info(x))
  ch.n.arg <- mappoly2:::embedded_to_numeric(chrom)
  assert_that(!any(is.na(ch.n.arg)), msg = "provide a valid chromosome identifier")
  ch.n.dat <- mappoly2:::embedded_to_numeric(x$chrom)
  ch.id <- which(ch.n.dat%in%ch.n.arg)
  if(has.mappoly2.screened(x))
    ch.id <- intersect(ch.id, get_screened_mrk_indices(x))
  ch.id
}

#' Given a dataset x and a vector containing the name of markers
#' w, returns metadata associated to w
get_sequence_metadata <- function(x,w){
  list(
    dosage.p1 = x$dosage.p1[w],
    dosage.p2 = x$dosage.p2[w],
    chrom = x$chrom[w],
    genome.pos = x$genome.pos[w],
    ref = x$ref[w],
    alt = x$alt[w],
    geno.dose = x$geno.dose[w],
    redundant = x$redundant[w],
    QAQC.values = x$QAQC.values$markers[w,]
  )
}

#' Removes all non-digit characters and then converts
#' the remaining string to an integer.
embedded_to_numeric <- function(x) {
  as.integer(gsub("[^0-9]", "", x))
}

#' Creates a formatted message with a horizontal rule
#' and specified text color.
#' @importFrom cli rule
#' @importFrom magrittr %>%
msg <- function(text, line = 1, col = "-"){
  cli::rule(line = line, left = text) %>%
    text_col(col = col) %>%
    message()
}

#' Changes the color of text in the R console, with a specific
#' focus on compatibility with RStudio's dark and light themes.
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




