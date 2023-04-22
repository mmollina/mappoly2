#' Pairwise Two-Point Analysis
#'
#' Performs a pairwise two-point analysis for all marker pairs in a given sequence.
#' The function estimates the recombination fraction for all possible linkage phase
#' configurations and their respective LOD Scores.
#'
#' @param input.seq An object of class \code{mappoly2.sequence}
#' @param ncpus Number of parallel processes (cores) to use (default = 1)
#' @param mrk.pairs A 2xN matrix containing N pairs of markers to analyze. If \code{NULL} (default), all pairs are considered
#' @param verbose If \code{TRUE} (default), displays information about the analysis
#' @param tol Desired accuracy. See \code{optimize()} for details
#' @param ll If \code{TRUE}, returns log-likelihood instead of LOD scores (for internal use)
#'
#' @return An object of class \code{mappoly2.twopt} which is a list containing the following components:
#' \item{input.seq}{The input sequence; an object of class \code{mappoly2.sequence} with the raw data}
#' \item{pairwise}{A list where each element is a data frame. Row names have the format x-y, where x and y represent the
#' number of homolos sharing the same allelic variant in parents P1 and P2, respectively (see Mollinari and Garcia, 2019 for notation).
#' The first column contains the LOD Score in relation to the most likely linkage phase configuration. The second column shows the estimated
#' recombination fraction for each configuration, and the third column contains the LOD Score comparing the likelihood under no linkage (r = 0.5)
#' with the estimated recombination fraction (evidence of linkage).}
#'
#' @export
#' @importFrom parallel makeCluster stopCluster parLapply
#' @importFrom utils combn
est_pairwise_rf <- function(input.seq,
                            ncpus = 1L,
                            mrk.pairs = NULL,
                            verbose = TRUE,
                            tol = .Machine$double.eps^0.25,
                            ll = FALSE)
{
  assert_that(is.mappoly2.sequence(input.seq))

  if (any(duplicated(input.seq$seq.num)))
    stop("There are duplicated markers in the sequence")
  if (is.null(mrk.pairs)) {
    seq.num <- sort(get_seq_indices(input.seq))
    mrk.pairs <- combn(sort(seq.num), 2) - 1
  } else {
    mrk.pairs <- mrk.pairs - 1
  }
  count.cache <- full_counts[[paste(input.seq$data[1:2], collapse = "x")]]
  ## splitting pairs in chunks
  if (length(input.seq$mrk.names) < 10)
    ncpus <- 1
  id <- ceiling(seq(1, (ncol(mrk.pairs) + 1), length.out = ncpus + 1))
  input.list <- vector("list", ncpus)
  for (i in 1:ncpus) input.list[[i]] <- mrk.pairs[, id[i]:(id[i + 1] - 1)]
  ## parallel version
  if (ncpus > 1) {
    start <- proc.time()
    if (verbose)
      cat("INFO: Using ", ncpus, " CPUs for calculation.\n")
    cl = parallel::makeCluster(ncpus)
    on.exit(parallel::stopCluster(cl))
    res <- parallel::parLapply(cl,
                               input.list,
                               paralell_pairwise,
                               input.seq = input.seq,
                               count.cache = count.cache,
                               tol = tol,
                               ll = ll)
    end <- proc.time()
    if (verbose) {
      cat("INFO: Done with",
          ncol(mrk.pairs),
          " pairs of markers \n")
      cat("INFO: Calculation took:",
          round((end - start)[3],
                digits = 3),
          "seconds\n")
    }
  }
  else {
    if (verbose) {
      cat("INFO: Going singlemode. Using one CPU for calculation.\n")
      if (length(input.seq$mrk.names) < 10)
        cat("Also, number of markers is too small to perform parallel computation.\n")
    }
    res <- lapply(input.list,
                  paralell_pairwise,
                  input.seq = input.seq,
                  count.cache = count.cache,
                  tol = tol,
                  ll = ll)
  }
  res <- unlist(res,
                recursive = FALSE)
  names(res) <- apply(mrk.pairs + 1,
                      2,
                      paste,
                      collapse = "-")
  nas <- sapply(res, function(x) any(is.na(x)))
  return(structure(list(input.seq = input.seq,
                        pairwise = res,
                        nas  = nas),
                   class = "mappoly2.twopt"))
}

print.mappoly2.twopt <- function(x)
  invisible(x)

#' Wrapper function to discrete-based pairwise two-point estimation in C++
#'
#' @param void internal function to be documented
#' @keywords internal
#' @export
paralell_pairwise <- function(mrk.pairs,
                              input.seq,
                              count.cache,
                              tol = .Machine$double.eps^0.25,
                              ll  = FALSE)
{
  res <- .Call("pairwise_rf_estimation",
               input.seq$data$ploidy.p1,
               input.seq$data$ploidy.p2,
               as.matrix(mrk.pairs),
               as.matrix(input.seq$data$geno.dose),
               as.vector(input.seq$data$dosage.p1),
               as.vector(input.seq$data$dosage.p2),
               count.cache,
               tol = tol,
               PACKAGE = "mappoly2")
  if(ll) return(res)
  return(lapply(res, format_rf))
  res
}

#' Format results from pairwise two-point estimation in C++
#'
#' @param void internal function to be documented
#' @keywords internal
format_rf <- function(res) {
  x <- res
  if (length(x) != 4) {
    LOD_ph <- min(x[2, ]) - x[2, ]
    rf <- x[1, order(LOD_ph, decreasing = TRUE)]
    LOD_rf <- abs(x[2, ] - x[3, ])[order(LOD_ph, decreasing = TRUE)]
    LOD_ph <- LOD_ph[order(LOD_ph, decreasing = TRUE)]
    return(cbind(LOD_ph, rf, LOD_rf))
  } else {
    nm <- paste(min(x[1], x[2]), min(x[3], x[4]), sep = "-")
    x <- cbind(0, rf = NA, LOD = NA)
    rownames(x) <- nm
    colnames(x) <- c("ph_LOD", "rf", "rf_LOD")
    return(x)
  }
}
