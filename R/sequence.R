
#' @export
make_sequence <- function(x,
                          lg = NULL,
                          ch = NULL,
                          mrk.id.list = NULL){
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
      return(mappoly2:::.map_skeleton(x$data, mrk.id.list))
    } else { # If not, use data
      if(!all(unlist(mrk.id.list) %in% x$screened.data$mrk.names))
        stop("Some markers names are not on the screened data")
      return(mappoly2:::.map_skeleton(x, mrk.id.list))
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
    return(mappoly2:::.map_skeleton(x$data, mrk.id.list))
  } else if(!is.null(lg) & is.null(ch)){
    assert_that(is.list(lg), msg = "'lg' should be a list")
    mrk.id.list <- vector("list", length(lg))
    for(i in 1:length(lg)){
      mrk.id.list[[i]] <- mappoly2:::get_markers_from_grouped_sequence(x, lg[[i]])
    }
    return(mappoly2:::.map_skeleton(x$data, mrk.id.list))
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
    return(mappoly2:::.map_skeleton(x$data, mrk.id.list))
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


