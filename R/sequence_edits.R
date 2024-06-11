#' Drop Marker(s) from a Sequence
#'
#' This function removes specified marker(s) from a given linkage group and type within a \code{mappoly2.sequence} object.
#' It can optionally re-estimate the map after marker removal.
#'
#' @param x A \code{mappoly2.sequence} object.
#' @param mrk Marker(s) to be dropped, specified by name or index.
#' @param lg Linkage group from which the marker(s) should be removed.
#' If \code{NULL} (default), all groups are considered.
#' @param type Type of sequence ('mds', 'genome', or 'custom') to specify which part of the sequence the marker(s) should be removed from.
#' @param parent Indicates the parent phase to consider ('p1p2', 'p1', or 'p2').
#' @param reestimate.map.and.haplo Logical; if \code{TRUE}, the genetic map and haplotype probabilities are re-estimated after marker removal (default is \code{FALSE}).
#' @param verbose Logical; if \code{TRUE}, function will print messages during execution (default is \code{TRUE}).
#' @param tol Tolerance level used in map re-estimation (applicable if \code{reestimate.map.and.haplo} is \code{TRUE}).
#'
#' @return Returns a modified \code{mappoly2.sequence} object with the specified marker(s) removed.
#'
#' @details The function identifies and removes the specified marker(s) from the linkage group and type within the \code{mappoly2.sequence} object.
#' If \code{reestimate.map.and.haplo} is \code{TRUE}, it also re-estimates the genetic map using the provided tolerance level.
#' The function validates the presence of the marker(s) and handles missing or incorrect marker specifications gracefully.
#'
#' @seealso \code{\link{add_marker}} for adding marker in a sequence.
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
drop_marker <- function(x,
                        mrk,
                        lg = NULL,
                        type = c("mds", "genome", "custom"),
                        parent = c("p1p2","p1","p2"),
                        reestimate.map.and.haplo = FALSE,
                        verbose = TRUE,
                        tol = 10e-4) {
  # Validate input object
  assert_that(is.mappoly2.sequence(x), msg = "Input 'x' must be a mappoly2.sequence object")

  # Ensure that 'mrk' can handle multiple markers
  if (!is.character(mrk)) {
    stop("Parameter 'mrk' must be a character vector of marker names")
  }

  y <- parse_lg_and_type(x, lg, type)
  parent <- match.arg(parent)

  # Check if the sequence is correctly mapped
  assert_that(is.mapped.sequence(x, y$lg, y$type, parent),
              msg = "The sequence is not correctly mapped for the given parameters")

  # Match markers and validate existence
  id <- match(mrk, rownames(x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$p1))
  if (any(is.na(id))) {
    if (verbose) {
      message("One or more specified markers are not in the original sequence. Returning original sequence.")
    }
    return(x)
  }

  # Confirm removal of markers if verbose
  if (verbose) {
    cat("Removing marker(s):", paste(mrk, collapse = ", "), "\n")
  }

  assert_that(is.logical(reestimate.map.and.haplo))
  mrk.names <- rownames(x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$p1)
  d <- cumsum(imf_h(c(0,x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$rf)))
  names(d) <- mrk.names
  if(any(!id%in%seq_along(mrk.names))){
    message("marker number(s) is not in within the range of the original sequence. \nreturning original sequence.")
    return(x)
  }
  x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$p1 <- x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$p1[-id, ]
  x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$p2 <-x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$p2[-id, ]
  x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$loglike <- NULL
  d <- mf_h(diff(d[-id]))
  names(d) <- NULL
  x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$rf <- d
  err <- x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$error
  x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$error <- NULL
  x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$haploprob <- NULL
  if(reestimate.map.and.haplo){
    x <- mapping(x, y$lg, y$type, parent, error = err,
                 rf = x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$rf,
                 tol = tol)
    x <- calc_haplotypes(x, y$lg, y$type, parent)
  }
  return(x)
}

#' Add a Marker to a Pre-Mapped Sequence
#'
#' This function adds a specified marker to a pre-mapped genetic sequence within a \code{mappoly2.sequence} object. It assesses the best position for the marker in the sequence based on phase configurations and recombination fractions.
#'
#' @param x A \code{mappoly2.sequence} object.
#' @param mrk The name of the marker to be added.
#' @param lg The linkage group to which the marker will be added.
#' @param type The type of sequence ('mds', 'genome', or 'custom').
#' @param parent Specifies the parent phase ('p1p2', 'p1', or 'p2').
#' @param reestimate.map.and.haplo Logical; if \code{TRUE}, re-estimates the map and haplotype probabilities after adding the marker (default is \code{FALSE}).
#' @param verbose Logical; if \code{TRUE}, function prints messages during execution (default is \code{TRUE}).
#' @param tol Tolerance level for re-estimating the map.
#' @param thresh.LOD.ph LOD threshold for considering phase configurations (default = 5).
#' @param thresh.LOD.rf LOD threshold for recombination fraction (default = 5).
#' @param thresh.rf Recombination fraction threshold (default = 0.5).
#' @param max.phases The maximum number of phase configurations to consider (default = 5).
#' @param thresh.LOD.ph.to.insert Threshold LOD for inserting the marker into the map (default = 10).
#' @param thresh.rf.to.add Recombination fraction threshold for adding the marker (default is the maximum of existing recombination fractions).
#'
#' @details
#' The function operates as follows:
#' 1. Validates the input object and parameters.
#' 2. Uses `test_one_marker` to evaluate all possible phase configurations for the specified marker.
#' 3. Checks if the best phase configuration meets the specified LOD and recombination fraction thresholds.
#' 4. Determines the optimal position for the marker within the linkage group based on the phase test results.
#' 5. Inserts the marker into the sequence at the determined position.
#' 6. Optionally re-estimates the map if `reestimate.map.and.haplo` is \code{TRUE}.
#'
#' @return Returns the modified \code{mappoly2.sequence} object with the new marker added.
#'
#' @seealso \code{\link{drop_marker}} for removing markers from a sequence.
#'
#' @export
add_marker <- function(x,
                       mrk,
                       lg = NULL,
                       type = c("mds", "genome", "custom"),
                       parent = c("p1p2","p1","p2"),
                       reestimate.map.and.haplo = FALSE,
                       verbose = TRUE,
                       tol = 10e-4,
                       thresh.LOD.ph = 5,
                       thresh.LOD.rf = 5,
                       thresh.rf = 0.5,
                       max.phases = 5,
                       thresh.LOD.ph.to.insert = 10,
                       thresh.rf.to.add = NULL){
  x1 <- x
  y <- parse_lg_and_type(x, lg, type)
  ph.res <- mappoly2:::test_one_marker(x,
                            mrk,
                            lg,
                            type,
                            parent,
                            thresh.LOD.ph,
                            thresh.LOD.rf,
                            thresh.rf,
                            max.phases,
                            verbose,
                            tol)
  if(length(ph.res) == 1)
    return(x1)

  if(length(ph.res$loglike) != 1 & diff(ph.res$loglike[1:2]) < thresh.LOD.ph.to.insert){
    if(verbose) message("LOD score below threshold.\nreturning original sequence.")
    return(x1)
  }
  if(is.null(thresh.rf.to.add))
    thresh.rf.to.add <- max(x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$rf)
  if(thresh.rf.to.add < 0 || thresh.rf.to.add >= 0.5)
    stop("'thresh.rf.to.add' argument must be between 0 and 0.5")
  if(max(ph.res$rf.vec[1,]) > thresh.rf.to.add){
    if(verbose)
      message("the recombination fraction excedes 'thresh.rf.to.add'\nreturning original seqeunce.")
    return(x1)
  }
  err <- x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$error
  pos <- ph.res$pos
  cur.mrk <- rownames(x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$p1)
  if(is.na(pos[[1]]$preceding))# beginning
  {
    x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$p1 <- rbind(ph.res$phases$p1[1,], x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$p1)
    x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$p2 <- rbind(ph.res$phases$p2[1,], x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$p2)
    rownames(x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$p1) <- rownames(x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$p2) <- c(mrk, cur.mrk)
    x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$rf <- c(ph.res$rf.vec[1,1], x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$rf)
  }
  else if (is.na(pos[[1]]$succeeding)){ #end
    x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$p1 <- rbind(x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$p1, ph.res$phases$p1[1,])
    x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$p2 <- rbind(x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$p2, ph.res$phases$p2[1,])
    rownames(x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$p1) <- rownames(x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$p2) <- c(cur.mrk, mrk)
    x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$rf <- c(x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$rf, ph.res$rf.vec[1,1])
  }
  else {
    idp <- 1:match(pos[[1]]$preceding, cur.mrk)
    ids <- (match(pos[[1]]$succeeding, cur.mrk)): length(cur.mrk)
    preceding <- cur.mrk[idp]
    succeeding <- cur.mrk[ids]
    x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$p1 <- rbind(x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$p1[preceding,],
                                                                  ph.res$phases$p1[1,],
                                                                  x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$p1[succeeding,])
    x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$p2 <- rbind(x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$p2[preceding,],
                                                                  ph.res$phases$p2[1,],
                                                                  x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$p2[succeeding,])
    rownames(x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$p1) <- rownames(x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$p2) <- c(preceding, mrk, succeeding)
    x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$rf <- c(x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$rf[idp[-length(idp)]],
                                                              ph.res$rf.vec[1,],
                                                              x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$rf[ids[-length(ids)]])
  }
  if(reestimate.map.and.haplo){
    x <- mapping(x, y$lg, y$type, parent, error = err,
                 rf = x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$rf,
                 tol = tol)
    x <- calc_haplotypes(x, y$lg, y$type, parent)
  }

  #s <- mappoly2:::calc_haploprob_biallelic_given_ve2()

  return(x)
}

#' Test All Possible Phases for a Given Marker in a Sequence
#'
#' This internal function tests all possible phase configurations for a specified marker within a linkage group of a \code{mappoly2.sequence} object. It is designed to evaluate how well a marker fits into an existing genetic map by exploring different phasing scenarios.
#'
#' @param x A \code{mappoly2.sequence} object.
#' @param mrk The name of the marker to be tested.
#' @param lg The linkage group to which the marker belongs.
#' @param type The type of sequence ('mds', 'genome', or 'custom').
#' @param parent Specifies the parent phase ('p1p2', 'p1', or 'p2').
#' @param thresh.LOD.ph LOD threshold for linkage phase (default = 5).
#' @param thresh.LOD.rf LOD threshold for recombination fraction (default = 5).
#' @param thresh.rf Recombination fraction threshold (default = 0.5).
#' @param max.phases The maximum number of phase configurations to consider (default = 5).
#' @param verbose Logical; if \code{TRUE}, function prints messages during execution (default is \code{TRUE}).
#' @param tol Tolerance level for phase estimation (default = 10e-4).
#'
#' @details
#' The function operates as follows:
#' 1. Validates the input object and the specified marker.
#' 2. Retrieves necessary information such as recombination fractions and dosage vectors for both parents.
#' 3. Performs two-point phasing for each parent to determine the most likely phase configurations.
#' 4. For each possible phase combination of the two parents, it estimates the Hidden Markov Model (HMM) map parameters.
#' 5. Compares different phasing scenarios based on the log-likelihood of the HMM model fit.
#' 6. The function returns a list of the best phase configurations based on the highest log-likelihood scores.
#'
#' @return A list containing the best phasing results, including log-likelihood scores, recombination fractions, and phase information for each parent.
#' @keywords internal
test_one_marker <- function(x,
                            mrk,
                            lg = 1,
                            type = c("mds", "genome", "custom"),
                            parent = c("p1p2","p1","p2"),
                            thresh.LOD.ph = 5,
                            thresh.LOD.rf = 5,
                            thresh.rf = 0.5,
                            max.phases = 5,
                            verbose = TRUE,
                            tol = 10e-4){
  assert_that(length(lg) == 1)

  # Validate input object
  assert_that(is.mappoly2.sequence(x), msg = "Input 'x' must be a mappoly2.sequence object")

  y <- mappoly2:::parse_lg_and_type(x, lg, type)
  parent <- match.arg(parent)

  # Check if the sequence is correctly mapped
  assert_that(is.haplotype.sequence(x, y$lg, y$type, parent),
              msg = "The sequence is not correctly mapped for the given parameters")
  assert_that(length(mrk)==1)
  assert_that(is.character(mrk))

  if(mrk%in%rownames(x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$p1)){
    message("'mrk' is alredy in the map. Retutning original sequence.")
    return(x)
  }
  ind.names <- x$data$screened.data$ind.names
  g <- x$data$geno.dose
  g[is.na(g)] <- -1
  mrk.pos <- rownames(x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$p1) # positioned markers
  mrk.seq <- mappoly2:::get_markers_from_ordered_sequence(x, y$lg, y$type, parent)[[1]] # all ordered marker
  assert_that(mrk %in% mrk.seq)
  assert_that(has.mappoly2.rf(x$data))

  M <- mappoly2:::filter_rf_matrix(x$data,
                        type = "sh",
                        thresh.LOD.ph,
                        thresh.LOD.rf,
                        thresh.rf,
                        mrk.names = mrk.seq)

  S1 <- M$Sh.p1[mrk, mrk.pos, drop = FALSE]
  S2 <- M$Sh.p2[mrk, mrk.pos, drop = FALSE]

  ## two-point phasing parent 1
  dose.vec1 <- x$data$dosage.p1[mrk]
  InitPh1 <- x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$p1
  L1 <- mappoly2:::phasing_one(mrk, dose.vec1, S1, InitPh1, FALSE)
  ## two-point phasing parent 2
  dose.vec2 <- x$data$dosage.p2[mrk]
  InitPh2 <- x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$p2
  L2 <- mappoly2:::phasing_one(mrk, dose.vec2, S2, InitPh2, FALSE)
  ## Selecting phase configurations
  n.conf <- sapply(L1, nrow) * sapply(L2, nrow)

  pedigree <- matrix(rep(c(1,
                           2,
                           x$data$ploidy.p1,
                           x$data$ploidy.p2, 1),
                         length(ind.names)),
                     nrow = length(ind.names),
                     byrow = TRUE)
  flanking <- mappoly2:::find_flanking_markers(mrk.seq, mrk.pos, mrk)
  G <- g[mrk, ind.names,drop = TRUE]
  u <- match(unlist(flanking[[mrk]]), mrk.pos)

  # Determine the position of the marker and set homolog probabilities
  if(is.na(u)[1]) {  # Marker at the beginning of the linkage group
    homolog_prob <- as.matrix(x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$haploprob[, c(na.omit(u), na.omit(u) + 1) + 3])
    idx <- c(0, 1, 2)
  } else if(is.na(u)[2]) {  # Marker at the end of the linkage group
    homolog_prob <- as.matrix(x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$haploprob[, c(na.omit(u) - 1, na.omit(u)) + 3])
    homolog_prob1 <- as.matrix(x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$haploprob[, na.omit(u) + 3, drop = FALSE])
    idx <- c(0, 2, 1)
  } else {  # Marker in the middle of the linkage group
    homolog_prob <- as.matrix(x$maps[[y$lg]][[y$type]][[parent]]$hmm.phase[[1]]$haploprob[, u + 3])
    idx <- c(1, 0, 2)
  }
  w2<-w1<-NULL
  z<-vector("list", nrow(L1[[1]]) * nrow(L2[[1]]))

  #if(length(z) == 1)
  #  phasing_results <- list(loglike = 0.000,
  #                          rf.vec = t(sapply(z[id],
  #                                            function(x) x[[2]])),
  #                          phases = list(p1 = w1[id,,drop=FALSE],
  #                                        p2 = w2[id,,drop=FALSE]),
  #                          pos = flanking)
  #phasing_results
  if(length(z) > max.phases)
    return(0)

  cat("~~~~~", length(z), "~~~~~\n")

  count <- 1
  for(j in 1:nrow(L1[[1]])){
    for(k in 1:nrow(L2[[1]])){
      PH <- list(L1[[1]][j,], L2[[1]][k,])
      if(is.na(u)[2]){
        z[[count]]<-mappoly2:::est_hmm_map_biallelic_insert_marker_at_the_end(PH,
                                                                    G,
                                                                    pedigree,
                                                                    homolog_prob1,
                                                                    rf = 0.01,
                                                                    verbose = FALSE,
                                                                    detailed_verbose = FALSE,
                                                                    tol = tol,
                                                                    ret_H0 = FALSE)
      } else {
        z[[count]]<-mappoly2:::est_hmm_map_biallelic_insert_marker(PH,
                                                                   G,
                                                                   pedigree,
                                                                   homolog_prob,
                                                                   rf = c(0.01,0.01),
                                                                   idx,
                                                                   verbose = FALSE,
                                                                   detailed_verbose = FALSE,
                                                                   tol = tol,
                                                                   ret_H0 = FALSE)
      }
      w1 <- rbind(w1, L1[[1]][j,])
      w2 <- rbind(w2, L2[[1]][k,])
      count <- count + 1
    }
  }
  v <- sapply(z, function(x) x[[1]])
  v <- max(v) - v
  id <- order(v)
  cat ("\n", v[id], "\n")
  phasing_results <- list(loglike = v[id],
                          rf.vec = t(sapply(z[id],
                                            function(x) x[[2]])),
                          phases = list(p1 = w1[id,,drop=FALSE],
                                        p2 = w2[id,,drop=FALSE]),
                          pos = flanking)
  phasing_results
}

#' Refine Genetic Map by Removing Problematic Markers
#'
#' Identifies and removes markers from a genetic map that contribute to significant gaps or inaccuracies,
#' potentially due to misplaced markers, incorrect phasing, or anomalous behavior.
#'
#' @param x A genetic mapping object, typically of a class "mappoly2.sequence".
#' @param lg Optional vector specifying the linkage group indices to process.
#'           If NULL, all linkage groups in `x` are processed.
#' @param type Character vector indicating the type of map to process, either "mds" or "genome".
#' @param parent Character vector specifying the parent or parents to be considered
#'               in the marker removal process. Options are "p1p2" (both parents),
#'               "p1" (first parent), and "p2" (second parent).
#' @param gap.threshold The threshold for identifying significant gaps in the map.
#'                      Markers causing gaps larger than this threshold will be considered for removal.
#' @param size.rem.cluster Minimum number of consecutive markers to retain.
#'                         Segments with fewer markers than this threshold will be removed.
#' @param ncpus The number of CPU cores to use for parallel processing. Defaults to 1.
#' @param reestimate.hmm.map If TRUE, re-estimates the hidden Markov model for the map after marker removal.
#' @param recompute.haplotype.prob If TRUE, recalculates haplotype probabilities after marker removal.
#' @param verbose Logical; if \code{TRUE}, function prints messages during execution (default is \code{TRUE}).
#' @param tol Tolerance level for the mapping algorithm.
#' @param error Optional; error rate to be used in the hidden Markov model estimation.
#'              If NULL, it uses the error rate from the input object.
#'
#' @return An updated genetic mapping object with problematic markers removed and, if specified,
#'         re-estimated HMM and recalculated haplotype probabilities.
#'
#' @details The function analyzes the genetic map to identify markers causing gaps
#'          larger than the specified threshold. It then removes these markers and, optionally,
#'          re-estimates the hidden Markov model and recalculates haplotype probabilities
#'          to refine the map.
#' @export
refine_map <- function(x,
                       lg = NULL,
                       type = c("mds", "genome"),
                       parent = c("p1p2", "p1", "p2"),
                       gap.threshold = 5,
                       size.rem.cluster =1,
                       ncpus = 1,
                       reestimate.hmm.map = TRUE,
                       recompute.haplotype.prob = TRUE,
                       verbose = TRUE,
                       tol = 10e-4,
                       error = NULL){
  # Extract the linkage group and type information from the input object
  y <- parse_lg_and_type(x, lg, type)

  assert_that(is.mappoly2.sequence(x))

  # Match the 'parent' argument to its possible values
  parent <- match.arg(parent)

  if(is.null(error))
    error <- x$maps[[1]][[y$type]][[parent]]$hmm.phase[[1]]$error

  cte <- logical(length(y$lg))
  for(i in y$lg){
    if(verbose)
      cat("Lg: ", i, "\n")
    assert_that(is.mapped.sequence(x, i, y$type, parent))
    adj.dist <- imf_h(x$maps[[i]][[y$type]][[parent]]$hmm.phase[[1]]$rf)
    n.mrk <- length(adj.dist) + 1
    id <- which( adj.dist > gap.threshold)
    if(length(id) == 0){
      cte[i] <- TRUE
      next()
    }
    id <- cbind(c(1, id+1), c(id, n.mrk))
    include.segments <- id[apply(id, 1, diff) > size.rem.cluster - 1, , drop = FALSE]
    remove.segments <-  id[apply(id, 1, diff) < size.rem.cluster, , drop = FALSE]
    if(nrow(remove.segments) == 0){
      cte[i] <- TRUE
      next()
    }
    if(length(include.segments) == 0) {
      cte[i] <- TRUE
      next()
    }
    rm <- NULL
    for(j in 1:nrow(remove.segments))
      rm <- c(rm, remove.segments[j,1]:remove.segments[j,2])
    mrk.rm <- rownames(x$maps[[i]][[y$type]][[parent]]$hmm.phase[[1]]$p1)[rm]
    x <- drop_marker(x, mrk = mrk.rm, lg = i, parent = parent,
                     type = y$type, verbose = verbose)
  }
  lg.temp <- lg[!cte]
  if(reestimate.hmm.map & any(!cte))
    x <- mapping(x,
                 lg = y$lg,
                 type = y$type,
                 parent = parent,
                 ncpus = ncpus,
                 error = error,
                 tol = tol)
  if(recompute.haplotype.prob & any(!cte))
    x <- calc_haplotypes(x,
                         lg = y$lg,
                         type = y$type,
                         parent = parent,
                         ncpus = ncpus)
  return(x)
}
