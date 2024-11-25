#' Create a Genetic Sequence from a Screened or Grouped Object
#'
#' This function generates a genetic sequence based on the provided screened or
#' grouped object. It allows for the creation of sequences using specific markers,
#' linkage groups, and chromosome information.
#'
#' @param x An object of class 'screened' or 'mappoly2.group'.
#' @param lg A list specifying the linkage groups, if applicable.
#' @param ch A list specifying the chromosomes, if applicable.
#' @param mrk.id.list A list or vector of marker IDs to be used in the sequence.
#'
#' @return Returns a sequence object based on the provided parameters.
#'         The sequence is constructed using internal function `.map_skeleton`.
#'
#' @details The function first checks if `mrk.id.list` is provided and is not null,
#'          and uses the marker information to initiate the sequence. If `lg` and `ch`
#'          are provided, they are used to determine the sequence structure.
#'          The function performs various checks and validations on the input parameters
#'          and the structure of the `x` object.
#'
#'          If `x` is of class 'mappoly2.group', the function uses data within the
#'          group object for sequence creation. If `lg` and `ch` are null, the function
#'          uses only the results of the UPGMA (Unweighted Pair Group Method with Arithmetic Mean).
#'
#' @importFrom assertthat assert_that
#' @export
make_sequence <- function(x,
                          lg = NULL,
                          ch = NULL,
                          mrk.id.list = NULL){

  assert_that(inherits(x, "mappoly2.data") | inherits(x, "mappoly2.group"))
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
      return(.map_skeleton(x$data, mrk.id.list))
    } else { # If not, use data
      if(!all(unlist(mrk.id.list) %in% x$screened.data$mrk.names))
        stop("Some markers names are not on the screened data")
      return(.map_skeleton(x, mrk.id.list))
    }
  }
  ## If lg and ch are null, use the results only of the UPGMA
  if(is.null(lg) & is.null(ch)){
    if(inherits(x, "mappoly2.group")){
      mrk.id.list <- split(names(x$groups.snp), x$groups.snp)
      return(.map_skeleton(x$data, mrk.id.list))
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
      mrk.id.list[[i]] <- get_markers_from_grouped_and_chromosome(x, lg[[i]], ch[[i]])
    }
    return(.map_skeleton(x$data, mrk.id.list))
  } else if(!is.null(lg) & is.null(ch)){
    assert_that(is.list(lg), msg = "'lg' should be a list")
    mrk.id.list <- vector("list", length(lg))
    for(i in 1:length(lg)){
      mrk.id.list[[i]] <- get_markers_from_grouped_sequence(x, lg[[i]])
    }
    return(.map_skeleton(x$data, mrk.id.list))
  } else if(is.null(lg) & !is.null(ch)){
    mrk.id.list <- vector("list", length(ch))
    assert_that(is.list(ch), msg = "'ch' should be a list")
    for(i in 1:length(ch)){
      mrk.id.list[[i]] <- x$data$mrk.names[get_mrk_indices_from_chrom(x$data, ch[[i]])]
      if(!is.null(x$initial.screened.rf))
        mrk.names[[i]] <- intersect(mrk.names[[i]], x$initial.screened.rf$mrk.names)
      else if(!is.null(x$screened.data$mrk.names)){
        mrk.id.list[[i]] <- intersect(mrk.names[[i]], x$screened.data$mrk.names)
      }
    }
    return(.map_skeleton(x$data, mrk.id.list))
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
    mrk.id.ch <- x$data$mrk.names[get_mrk_indices_from_chrom(x$data, ch)]
  }
  assert_that(length(mrk.id.ch) > 0)
  if(!is.null(lg))
    mrk.id.lg <- get_markers_from_grouped_sequence(x, lg)
  return(intersect(mrk.id.ch, mrk.id.lg))
}


