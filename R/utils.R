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
  a <- unique(x$chrom)
  na.name <- a[is.na(mappoly2:::embedded_to_numeric(a))]
  ch.n.arg <- mappoly2:::embedded_to_numeric(chrom)
  if(any(is.na(ch.n.arg))){
    b <- chrom%in%na.name & is.na(ch.n.arg)
    if(all(!b))
      stop(msg = "provide a valid chromosome identifier")
  }
  ch.n.dat <- mappoly2:::embedded_to_numeric(x$chrom)
  ch.n.arg[is.na(ch.n.arg)] <- ch.n.dat[is.na(ch.n.dat)] <- 999
  ch.id <- which(ch.n.dat%in%ch.n.arg)
  if(has.mappoly2.screened(x))
    ch.id <- intersect(ch.id, get_screened_mrk_indices(x))
  ch.id
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

#' Sort Phases Based on Log-Likelihood
#'
#' This function sorts a list of phases based on their log-likelihood values.
#' It can return either the best phase (highest log-likelihood) or all phases in a sorted order.
#'
#' @param ph A list of phases, each phase should have a `$loglike` element.
#' @param only.best A boolean flag; if TRUE, only the best (highest log-likelihood) phase is returned.
#'                   If FALSE, all phases are returned in descending order of their log-likelihoods. Default is TRUE.
#' @return If `only.best` is TRUE, returns the phase with the highest log-likelihood.
#'         If FALSE, returns all phases sorted by log-likelihood in descending order.
sort_phase <- function(ph, only.best = TRUE){
  ph <- ph[sapply(ph, function(x) !is.null(x$loglike))]
  if(only.best){
    ph <- ph[which.max(sapply(ph, function(x) x$loglike))]
    return(ph)
  } else {
    ph <- ph[order(sapply(ph, function(x) x$loglike), decreasing = TRUE)]
    return(ph)
  }
}

#' Find Flanking Markers
#'
#' Identifies the closest preceding and succeeding elements from `B1` for each element in `B2`,
#' based on their positions in `A`.
#'
#' @param A A character vector; positions of elements in this vector are considered.
#' @param B1 A character vector; contains reference elements.
#' @param B2 A character vector; elements for which closest flanking elements in `B1` are to be found.
#' @return A list where each element in `B2` is mapped to its closest preceding and succeeding elements from `B1`.
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



#' Parse Linkage Group and Type
#'
#' This function parses the linkage group and type from a given mappoly2.sequence object.
#' It validates the input and extracts the specified linkage group and type.
#'
#' @param x An object of class mappoly2.sequence. This object should contain SNP data
#'   along with associated maps.
#' @param lg A character or numeric vector specifying the linkage groups. If a character
#'   vector, it should correspond to the names of the maps in the mappoly2.sequence object.
#'   If numeric, it should be the indices of the maps.
#' @param type A character vector specifying the type of SNP order.
#'   Possible values are "mds", "genome", or "custom". Default is c("mds", "genome", "custom").
#'   The function will match the provided type argument with these options.
#'
#' @return A list containing two elements: 'lg' and 'type'. The 'lg' element contains the
#'   linkage group information, either as a character or numeric vector, as provided by the user.
#'   The 'type' element contains the matched SNP order type.
#'
#' @export
#'
#' @importFrom assertthat assert_that
parse_lg_and_type <- function(x, lg = NULL, type = c("mds", "genome", "custom")){
  assert_that(is.mappoly2.sequence(x))
  type <- match.arg(type)
  if(is.null(lg))
    lg <- seq_along(x$maps)
  if(is.character(lg)){
    assert_that(lg %in% names(x$maps), msg = "Provide a valid group set")
    lg <- match(lg, names(x$maps))
  }
  assert_that(all(lg <= length(x$maps)))
  return(list(lg = lg, type = type))
}


#' Extract Markers from Ordered Sequence
#'
#' Extracts markers from an ordered sequence in an object, based on the specified type.
#'
#' @param x An object containing marker maps.
#' @param lg Either a character vector naming the maps or a numeric vector indicating map indices.
#' @param type The type of sequence to consider (either "mds", "genome", or "custom").
#' @return Markers based on the specified type.
get_markers_from_ordered_sequence <- function(x, lg, type = c("mds", "genome", "custom")){
  y <- mappoly2:::parse_lg_and_type(x,lg,type)
  if(y$type == "mds")
    return(lapply(x$maps[y$lg], function(z, type) intersect(z[[type]]$order$locimap$locus,
                                                   z[[type]]$mkr.names), y$type))
  else if(y$type == "genome")
    return(lapply(x$maps[y$lg], function(z, type) intersect(rownames(z[[type]]$order),
                                                     z[[type]]$mkr.names), y$type))
}

#' Extract Markers from Mapped Sequence
#'
#' Similar to `get_markers_from_ordered_sequence` but specifically for mapped sequences.
#'
#' @param x An object containing marker maps.
#' @param lg Either a character vector naming the maps or a numeric vector indicating map indices.
#' @param type The type of sequence to consider (either "mds", "genome", or "work").
#' @return Markers based on the specified type for mapped sequences.
get_markers_from_phased_sequence <- function(x, lg, type = c("mds", "genome", "custom")){
  y <- mappoly2:::parse_lg_and_type(x,lg,type)
  lapply(x$maps[y$lg], function(z, type) intersect(rownames(z[[type]]$phase[[1]]$p1),
                                                   z[[type]]$mkr.names), y$type)
}

#' Determine Dosage Type of Markers
#'
#' Classifies each marker based on its dosage type using the provided dosage data.
#'
#' @param x An object containing dosage data.
#' @param mrk.names Names of the markers to be analyzed.
#' @return A list categorizing each marker as simplex in parent 1, simplex in parent 2,
#'         double simplex, or multiplex based on their dosage values.
get_dosage_type <- function(x, mrk.names){
  p1 <- abs(abs(x$data$dosage.p1 - x$data$ploidy.p1/2) - x$data$ploidy.p1/2)
  p2 <- abs(abs(x$data$dosage.p2 - x$data$ploidy.p2/2) - x$data$ploidy.p2/2)
  s.p1 <- p1  ==  1 & p2  ==  0
  s.p2 <- p1  ==  0 & p2  ==  1
  ds <- p1  ==  1 & p2  ==  1
  list(simplex.p1 = names(which(s.p1[mrk.names])),
       simplex.p2 = names(which(s.p2[mrk.names])),
       double.simplex = names(which(ds[mrk.names])),
       multiplex = names(which(!(s.p1 | s.p2 | ds)[mrk.names])))
}

#' Create a Map Skeleton for mappoly2 Data
#'
#' This function creates a skeleton structure for mapping data, specifically for mappoly2 data objects.
#' It initializes map structures for multiple linkage groups based on provided marker IDs.
#'
#' @param x A mappoly2 data object; the function checks if the input is a valid mappoly2 data structure.
#' @param mrk.id.list A list or a vector of marker IDs for which the map skeleton is to be created.
#'                    If a vector is provided, it is converted into a list.
#' @return Returns a mappoly2 sequence object with initialized map structures for each linkage group.
#'         Each linkage group will have a map for MDS, genome, and work, but with `NULL` values for `order` and `phase`.
#' @note This function is internal and usually not called directly by the user.
.map_skeleton<- function(x, mrk.id.list){
  assert_that(is.mappoly2.data(x))
  if(!is.list(mrk.id.list))
    mrk.id.list <- list(mrk.id.list)
  n.lg <- length(mrk.id.list)
  maps <- vector("list", n.lg)
  names(maps) <- c(paste0("lg", 1:n.lg))
  for(i in 1:n.lg){
    maps[[i]] <- list(mds = list(mkr.names = mrk.id.list[[i]],
                                 order = NULL,
                                 phase = NULL),
                      genome = list(mkr.names = mrk.id.list[[i]],
                                    order = NULL,
                                    phase = NULL),
                      custom = list(mkr.names = mrk.id.list[[i]],
                                    order = NULL,
                                    phase = NULL))
  }
  structure(list(maps = maps, data = x), class = "mappoly2.sequence")
}
