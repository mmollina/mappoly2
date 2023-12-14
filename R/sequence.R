#' Create a sequence of markers
#'
#' Makes a sequence of markers based on an object of another class.
#'
#' @param x an object of one of the following classes:
#'     \code{mappoly.data}, \code{mappoly.map}, \code{mappoly.sequence},
#'     \code{mappoly.group},
#'     \code{mappoly.pcmap}, \code{mappoly.pcmap3d}, or \code{mappoly.geno.ord}
#'
#' @param arg one of the following input types:
#' \enumerate{
#' \item A string 'all': Generates a sequence with all markers from the raw data.
#' \item A string or a vector of strings 'chrX': Specifies a chromosome (with 'X'
#'       being the chromosome number). For unassigned markers, use 'chr0'.
#' \item A vector of integers: Indicates the position of markers in the original
#'       data set to be included in the sequence.
#' \item A vector of strings: Indicates the names or identifiers of the genetic
#'       markers
#' \item An integer: Represents a linkage group when xect is of class
#'       \code{mappoly.group}.
#' \item NULL: Applicable when \code{xect} belongs to one of the following
#'       classes: \code{mappoly.pcmap}, \code{mappoly.pcmap3d},
#'       \code{mappoly.unique.seq}, or \code{mappoly.geno.ord}.
#' }
#'
#' @param info.parent one of the following options:
#' \enumerate{
#' \item \code{'all'}{select all dosage combinations in both parents (default)}
#' \item \code{'P1'}{select informative markers parent 1}
#' \item \code{'P2'}{select informative markers parent 2}
#' }
#'
#' @param genomic.info An optional argument applicable only to \code{mappoly.group}
#' objects, which can be either NULL or a numeric combination of sequences
#' from genomic information used to create the sequences:
#' \enumerate{
#' \item NULL (default): Returns a sequence containing all markers as defined by
#' the grouping function.
#' \item 1: Returns a sequence with markers that match the intersection between
#' the grouping function and genomic information, considering the genomic information
#' sequence with the maximum number of matching markers for the group.
#' \item c(1, 2): Returns a sequence with markers that match the intersection between
#' the grouping function and genomic information, considering the two genomic
#' information sequences with the maximum number of matching markers for the group,
#' and so on.
#' }
#'
#' @param x an object of the class \code{mappoly.sequence}
#'
#' @param thresh.line position of a threshold line for p values of the segregation test (default = \code{0.05/n.mrk})
#'
#' @param ... currently ignored
#'
#' @return An object of class \code{mappoly2.sequence}, which is a
#'     list containing the following components:
#'     \item{mrk.names}{}
#'     \item{phases}{}
#'     \item{redundant}{}
#'     \item{data}{}
#'
#' @author Marcelo Mollinari (\email{mmollin@ncsu.edu}) and
#'         Gabriel Gesteira, (\email{gdesiqu@ncsu.edu})
#' @export
#' @importFrom assertthat assert_that
#'
# make_sequence <- function(input.data,
#                           grouping.scope = c("lg.grouping",
#                                              "chrom",
#                                              "mrk.seq"),
#                           genomic.info = NULL)
#   assert_that(has.mappoly2.rf(input.data))
#
# #### Grouping ####
# x <- group(x, expected.groups = 13)
# plot(x, what = "group")
# print(x)
# #### Set Working Sequences ####
# x <- set_working_sequence(x,
#                           lg = list(c(1,4),
#                                     c(2,3),
#                                     c(5,6)),
#                           ch = list(c(1,3),
#                                     2,
#                                     c(4,5)))
make_sequence <- function(x,
                          lg = NULL,
                          ch = NULL,
                          mrk.id.list = NULL,
                          ind.names = NULL,
                          seq.names = NULL,
                          seq.order = NULL){
  assert_that(inherits(x, "screened") | inherits(g, "mappoly2.group"))
  ## If there is marker information, initiate a sequence with it
  ## It over hides lg and ch
  if(!is.null(mrk.id.list)){
    if(!is.list(mrk.id.list)){
      mrk.id.list <- list(mrk.id.list)
    }
    if(any(!sapply(mrk.id.list, is.character)))
      stop("Provide the name of the markers, not their indices")
    if(inherits(x, "mappoly2.group")){ # If grouped, use the data within group object
      if(!all(unlist(mrk.id.list) %in% x$data$screened.data$mrk.names))
        stop("Some markers names are not on the screened data")
      if(is.null(seq.order)) seq.order <- 1:length(mrk.id.list)
      return(mappoly2:::.map_skeleton(x$data, mrk.id.list[seq.order]))
    } else { # If not, use data
      if(!all(unlist(mrk.id.list) %in% x$screened.data$mrk.names))
        stop("Some markers names are not on the screened data")
      if(is.null(seq.order)) seq.order <- 1:length(mrk.id.list)
      return(mappoly2:::.map_skeleton(x, mrk.id.list[seq.order]))
    }
  }
  ## If lg and ch are null, use the results only of the UPGMA
  if(is.null(lg) & is.null(ch)){
    if(inherits(x, "mappoly2.group")){
      mrk.id.list <- split(names(x$groups.snp), x$groups.snp)
      return(mappoly2:::.map_skeleton(x$data, mrk.id.list))
    } else {
      stop("Provide a list of markers of a grouped object")
    }
  }
  else if(!is.null(lg) & !is.null(ch)){
    assert_that(is.list(lg), msg = "'lg' should be a list")
    assert_that(is.list(ch), msg = "'ch' should be a list")
    assert_that(length(lg) == length(ch))
    mrk.id.list <- vector("list", length(lg))
    for(i in 1:length(lg)){
      mrk.id.list[[i]] <- mappoly2:::get_markers_from_grouped_and_chromosome(x, lg[[i]], ch[[i]])
    }
    if(is.null(seq.order)) seq.order <- 1:length(mrk.id.list)
    return(mappoly2:::.map_skeleton(x$data, mrk.id.list[seq.order]))
  } else if(!is.null(lg) & is.null(ch)){
    assert_that(is.list(lg), msg = "'lg' should be a list")
    mrk.id.list <- vector("list", length(lg))
    for(i in 1:length(lg)){
      mrk.id.list[[i]] <- mappoly2:::get_markers_from_grouped_sequence(x, lg[[i]])
    }
    if(is.null(seq.order)) seq.order <- 1:length(mrk.id.list)
    return(mappoly2:::.map_skeleton(x$data, mrk.id.list[seq.order]))
  } else if(is.null(lg) & !is.null(ch)){
    mrk.id.list <- vector("list", length(ch))
    assert_that(is.list(ch), msg = "'ch' should be a list")
    for(i in 1:length(ch)){
      mrk.id.list[[i]] <- x$data$mrk.names[mappoly2:::get_mrk_indices_from_chrom(x$data, ch[[i]])]
      if(!is.null(x$initial.screened.rf))
        mrk.names[[i]] <- intersect(mrk.names[[i]], x$initial.screened.rf$mrk.names)
      else if(!is.null(x$screened.data$mrk.names)){
        mrk.id.list[[i]] <- intersect(mrk.names[[i]], x$screened.data$mrk.names)
      }
    }
    if(is.null(seq.order)) seq.order <- 1:length(mrk.id.list)
    return(mappoly2:::.map_skeleton(x$data, mrk.id.list[seq.order]))
  }
}

get_markers_from_grouped_sequence <- function(x, arg){
  inherits(x, "mappoly2.group")
  assert_that(is.numeric(arg))
  names(x$groups.snp)[x$groups.snp  %in%  arg]
}

get_markers_from_grouped_and_chromosome <- function(x, lg = NULL, ch = NULL){
  inherits(x, "mappoly2.group")
  assert_that(!is.null(lg) | !is.null(ch),
              msg = "Please provide a value for either 'lg' or 'ch'. Both cannot be left blank.")
  mrk.id.ch <- NULL
  if(!is.null(ch)){
    mrk.id.ch <- x$data$mrk.names[mappoly2:::get_mrk_indices_from_chrom(x$data, ch)]
  }
  assert_that(length(mrk.id.ch) > 0)
  if(!is.null(lg))
    mrk.id.lg <- mappoly2:::get_markers_from_grouped_sequence(x, lg)
  return(intersect(mrk.id.ch, mrk.id.lg))
}


