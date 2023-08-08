#' @export pairwise_rf
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
pairwise_rf <- function(x,
                        ncpus = 1L,
                        mrk.pairs = NULL,
                        verbose = TRUE,
                        tol = .Machine$double.eps^0.25)
{
 assert_that(inherits(x, "screened"),
             msg = "Error: The input object has not been screened.")
  if (is.null(mrk.pairs)) {
    seq.num <- seq_along(x$screened.data$mrk.names)
    mrk.pairs <- combn(sort(seq.num), 2) - 1
  } else {
    mrk.pairs <- mrk.pairs - 1
  }
  count.cache <- mappoly2:::full_counts[[paste(sort(unlist(x$data[1:2])), collapse = "x")]]
  RcppParallel::setThreadOptions(numThreads = ncpus)
  count.vector = unlist(count.cache)
  count.phases = unlist(lapply(count.cache, function(x) paste0(names(x), collapse = '/')))
  count.matrix.rownames = unlist(lapply(count.cache, function(x) paste0(rownames(x[[1]]), collapse = '/')))
  count.matrix.number = unlist(lapply(count.cache, length))
  count.matrix.length = unlist(lapply(count.cache, function(x) length(c(unlist(x)) )))
  count.matrix.pos = cumsum(c(1, count.matrix.length[-length(count.matrix.length)]))
  ploidy.p1 <- x$data$ploidy.p1
  ploidy.p2 <- x$data$ploidy.p2
  if(ploidy.p1 <= ploidy.p2){
    dose.p1 <- x$data$dosage.p1[x$screened.data$mrk.names]
    dose.p2 <- x$data$dosage.p2[x$screened.data$mrk.names]
    swap.parents <- FALSE
  } else {
    ploidy.p1 <- x$data$ploidy.p2
    ploidy.p2 <- x$data$ploidy.p1
    dose.p1 <- x$data$dosage.p2[x$screened.data$mrk.names]
    dose.p2 <- x$data$dosage.p1[x$screened.data$mrk.names]
    swap.parents <- TRUE
  }
  geno <- as.matrix(x$data$geno.dose[x$screened.data$mrk.names, x$screened.data$ind.names])
  geno[is.na(geno)] <- 1 + (ploidy.p1 + ploidy.p2)/2
  res <- mappoly2:::pairwise_rf_estimation_disc_rcpp(mrk_pairs_R = as.matrix(mrk.pairs),
                                          ploidy_p1_R = ploidy.p1,
                                          ploidy_p2_R = ploidy.p2,
                                          geno_R = geno,
                                          dose_p1_R = as.vector(dose.p1),
                                          dose_p2_R = as.vector(dose.p2),
                                          count_vector_R = count.vector,
                                          count_matrix_phases_R = count.phases,
                                          count_matrix_rownames_R = count.matrix.rownames,
                                          count_matrix_number_R = count.matrix.number,
                                          count_matrix_pos_R = count.matrix.pos,
                                          count_matrix_length_R = count.matrix.length,
                                          tol_R = tol, threads_R = ncpus)
  res[res == -1] = NA
  colnames(res) = c("Sh_P1","Sh_P2","rf","LOD_rf","LOD_ph")
  seq.mrk.names <- x$screened.data$mrk.names
  v_2_m <- function(x, n){
    y <- base::matrix(NA, n, n)
    y[base::lower.tri(y)] <- as.numeric(x)
    y[base::upper.tri(y)] <- t(y)[base::upper.tri(y)]
    y
  }
  n <- length(seq.mrk.names)
  if(swap.parents){
    Sh_P2 = v_2_m(res[,1], n)
    Sh_P1 = v_2_m(res[,2], n)
  } else {
    Sh_P1 = v_2_m(res[,1], n)
    Sh_P2 = v_2_m(res[,2], n)
  }
  rf = v_2_m(res[,3], n)
  LOD_rf = v_2_m(res[,4], n)
  LOD_ph = v_2_m(res[,5], n)
  dimnames(rf) <- dimnames(LOD_rf) <- dimnames(LOD_ph) <- dimnames(Sh_P1) <- dimnames(Sh_P2) <- list(seq.mrk.names, seq.mrk.names)
  x$pairwise = list(rec.mat = rf,
                            lod.mat = abs(LOD_rf),
                            lod.ph.mat = abs(LOD_ph),
                            Sh.p1 = Sh_P1,
                            Sh.p2 = Sh_P2)
  class(x) <- c(class(x), "pairwise")
  return(x)
}



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
#' number of homologs sharing the same allelic variant in parents P1 and P2, respectively (see Mollinari and Garcia, 2019 for notation).
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
  count.cache <- full_counts[[paste(sort(unlist(input.seq$data[1:2])), collapse = "x")]]
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
  ploidy.p1 <- input.seq$data$ploidy.p1
  ploidy.p2 <- input.seq$data$ploidy.p2
  if(ploidy.p1 <= ploidy.p2){
    dose.p1 <- input.seq$data$dosage.p1
    dose.p2 <- input.seq$data$dosage.p2
    swap.parents <- FALSE
  } else {
    ploidy.p1 <- input.seq$data$ploidy.p2
    ploidy.p2 <- input.seq$data$ploidy.p1
    dose.p1 <- input.seq$data$dosage.p2
    dose.p2 <- input.seq$data$dosage.p1
    swap.parents <- TRUE
  }
  res <- .Call("pairwise_rf_estimation",
               ploidy.p1,
               ploidy.p2,
               as.matrix(mrk.pairs),
               as.matrix(input.seq$data$geno.dose),
               as.vector(dose.p1),
               as.vector(dose.p2),
               count.cache,
               tol,
               swap.parents,
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

#' Recombination fraction list to matrix
#'
#' Transforms the recombination fraction list contained in an object
#' of class \code{mappoly.twopt} or \code{mappoly.twopt2} into a recombination
#' fraction matrix
#'
#' \code{thresh_LOD_ph} should be set in order to only select
#'     recombination fractions that have LOD scores associated to the
#'     linkage phase configuration higher than \code{thresh_LOD_ph}
#'     when compared to the second most likely linkage phase configuration.
#'
#' @param input.twopt an object of class \code{mappoly.twopt} or \code{mappoly.twopt2}
#'
#' @param thresh.LOD.ph LOD score threshold for linkage phase configurations (default = 0)
#'
#' @param thresh.LOD.rf LOD score threshold for recombination fractions (default = 0)
#'
#' @param thresh.rf the threshold used for recombination fraction filtering (default = 0.5)
#'
#' @param ncpus number of parallel processes (i.e. cores) to spawn (default = 1)
#'
#' @param shared.alleles if \code{TRUE}, computes two matrices (for both parents) indicating
#'                       the number of homologues that share alleles (default = FALSE)
#'
#' @param verbose if \code{TRUE} (default), current progress is shown; if
#'     \code{FALSE}, no output is produced
#'
#' @param x an object of class \code{mappoly.rf.matrix}
#'
#' @param type type of matrix that should be printed. Can be one of the
#'        following: \code{"rf"}, for recombination fraction or \code{"lod"}
#'        for LOD Score
#'
#' @param ord the order in which the markers should be plotted (default = NULL)
#'
#' @param rem which markers should be removed from the heatmap (default = NULL)
#'
#' @param main.text a character string as the title of the heatmap (default = NULL)
#'
#' @param index \code{logical} should the name of the markers be printed in the
#' diagonal of the heatmap? (default = FALSE)
#'
#' @param fact positive integer. factor expressed as number of cells to be aggregated
#' (default = 1, no aggregation)
#'
#' @param ... currently ignored
#'
#' @return A list containing two matrices. The first one contains the
#'     filtered recombination fraction and the second one contains the
#'     information matrix
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_.
#'     \doi{10.1534/g3.119.400378}
#'
#' @importFrom fields tim.colors image.plot
#' @export rf_list_to_matrix

rf_list_to_matrix <- function(input.twopt,
                              thresh.LOD.ph = 0,
                              thresh.LOD.rf = 0,
                              thresh.rf = 0.5,
                              ncpus = 1L,
                              shared.alleles = FALSE,
                              verbose = TRUE) {
  ## checking for correct object
  assert_that(is.mappoly2.twopt(input.twopt))
  pair_input <- input.twopt$pairwise
  n.mrk <- length(input.twopt$input.seq$mrk.names)
  seq.num <- sort(get_seq_indices(input.twopt$input.seq))
  marnames <- rownames(input.twopt$input.seq$data$geno.dose)[seq.num]
  lod.ph.mat <- lod.mat <- rec.mat <- matrix(NA, n.mrk, n.mrk)
  if(shared.alleles)
    Sh.p1 <- Sh.p2 <- lod.mat

  #### UPDATE: instead of recovering the order from names, provide using the object 'input.twopt'
  #seq.num.orig <- unique(sapply(strsplit(x = names(input.twopt$pairwise), split = "-"), function(x) as.numeric(x[1])))
  #seq.num.orig <- c(seq.num.orig, strsplit(tail(names(input.twopt$pairwise), n = 1), "-")[[1]][2])
  #dimnames(lod.mat) = dimnames(rec.mat) = list(seq.num.orig, seq.num.orig)
  if (ncpus > 1) {
    start <- proc.time()
    if (verbose)
      cat("INFO: Using ", ncpus, " CPUs.\n")
    cl <- parallel::makeCluster(ncpus)
    parallel::clusterExport(cl,
                            varlist = c("thresh.LOD.ph", "thresh.LOD.rf", "thresh.rf", "shared.alleles"),
                            envir = environment())
    rf.lod.mat <- parallel::parSapply(cl, pair_input, "select_rf",
                                      thresh.LOD.ph, thresh.LOD.rf, thresh.rf, shared.alleles)
    parallel::stopCluster(cl)
    end <- proc.time()
    if (verbose) {
      cat("INFO: Done with",
          length(pair_input),
          " pairs of markers \n")
      cat("INFO: Operation took:",
          round((end - start)[3],
                digits = 3),
          "seconds\n")
    }
  } else {
    if (verbose) {
      cat("INFO: Going singlemode. Using one CPU.\n")
    }
    rf.lod.mat <- sapply(pair_input, select_rf, thresh.LOD.ph, thresh.LOD.rf, thresh.rf, shared.alleles)
  }
  rec.mat[lower.tri(rec.mat)] <- as.numeric(rf.lod.mat[1,])
  rec.mat[upper.tri(rec.mat)] <- t(rec.mat)[upper.tri(rec.mat)]
  lod.mat[lower.tri(lod.mat)] <- as.numeric(rf.lod.mat[2,])
  lod.mat[upper.tri(lod.mat)] <- t(lod.mat)[upper.tri(lod.mat)]
  lod.ph.mat[lower.tri(lod.ph.mat)] <- as.numeric(rf.lod.mat[3,])
  lod.ph.mat[upper.tri(lod.ph.mat)] <- t(lod.ph.mat)[upper.tri(lod.ph.mat)]
  dimnames(lod.ph.mat) <- dimnames(rec.mat) <- dimnames(lod.mat) <- list(marnames, marnames)
  if(shared.alleles){
    Sh.p1[lower.tri(Sh.p1)] <- as.numeric(rf.lod.mat[4,])
    Sh.p1[upper.tri(Sh.p1)] <- t(Sh.p1)[upper.tri(Sh.p1)]
    Sh.p2[lower.tri(Sh.p2)] <- as.numeric(rf.lod.mat[5,])
    Sh.p2[upper.tri(Sh.p2)] <- t(Sh.p2)[upper.tri(Sh.p2)]
    dimnames(Sh.p1) <- dimnames(Sh.p2) <- list(marnames, marnames)
  } else{
    Sh.p1 <- Sh.p2 <- NULL
  }
  structure(list(thresh.LOD.ph = thresh.LOD.ph,
                 thresh.LOD.rf = thresh.LOD.rf,
                 thresh.rf = thresh.rf,
                 rec.mat = rec.mat,
                 lod.mat = abs(lod.mat),
                 lod.ph.mat = abs(lod.ph.mat),
                 Sh.p1 = Sh.p1,
                 Sh.p2 = Sh.p2,
                 input.seq = input.twopt$input.seq),
            class = "mappoly2.rf.matrix")
}

#' @rdname rf_list_to_matrix
#' @export
plot_mappoly2_rf_matrix <- function(x, type = c("rf", "lod"), ord = NULL, rem = NULL,
                                    main.text = NULL, index = FALSE, fact = 1, ...){
  type <- match.arg(type)
  if(is.mappoly2.sequence(ord))
    ord <- ord$mrk.names
  if(type  ==  "rf"){
    w <- x$rec.mat
    if(!is.null(ord))
    {
      w <- w[ord,ord]
    }
    if(!(is.null(rem) || sum(colnames(x$rec.mat)%in%rem)  ==  0))
    {
      o <- which(colnames(x$rec.mat)%in%rem)
      w <- w[-o,-o]
    }
    if(fact > 1)
      w <- aggregate_matrix(w, fact)
    if(is.null(main.text))
      main.text <- "Recombination fraction matrix"
    col.range  <-
      na.omit(rev(fields::tim.colors())[1:(ceiling(128 * max(x$rec.mat, na.rm = TRUE)) + 1)])
    brks <- NULL
  } else if(type  ==  "lod")
  {
    w <- x$lod.mat
    if(!is.null(ord))
    {
      w <- w[ord,ord]
    }
    if(!(is.null(rem) || sum(colnames(x$rec.mat)%in%rem)  ==  0))
    {
      o <- which(colnames(x$rec.mat)%in%rem)
      w <- w[-o,-o]
    }
    if(fact > 1)
      w <- aggregate_matrix(w, fact)
    w[w < 1e-4] <- 1e-4
    w <- log10(w)
    if(is.null(main.text))
      main.text <- "log(LOD) Score matrix"
    col.range <- na.omit(fields::tim.colors()[1:(ceiling(128 * max(x$lod.mat, na.rm = TRUE)) + 1)])
    col.range <- col.range[ceiling(seq(1, length(col.range), length.out = 10))]
    brks <- seq(min(w, na.rm = TRUE), max(w, na.rm = TRUE), length.out = 11)
    brks <- round(exp(brks/log10(exp(1))),1)
  } else stop("Invalid matrix type.")

  fields::image.plot(
    w,
    col = col.range,
    lab.breaks = brks,
    main = main.text,
    useRaster = FALSE,
    axes = FALSE
  )
  if(ncol(w) < 100)
    ft <- .7
  else
    ft <- 100/ncol(w)
  if(index)
    text(x = seq(0,1, length.out = ncol(w)), y = seq(0,1, length.out = ncol(w)),
         labels = colnames(w), cex = ft)
}


#' Select rf and lod based on thresholds
#'
#' @param void internal function to be documented
#' @keywords internal
select_rf <- function(x, thresh.LOD.ph, thresh.LOD.rf, thresh.rf, shared.alleles = FALSE)
{
  if(any(is.na(x)))
  {
    if(shared.alleles){return(c(NA,NA,NA,NA,NA))} else return(c(NA,NA,NA))
  }
  if((nrow(x)  ==  1 || abs(x[2 , 1]) >= thresh.LOD.ph) &&
     abs(x[1, 3]) >= thresh.LOD.rf &&
     abs(x[1, 2]) <= thresh.rf)
  {
    if(shared.alleles){
      y <- strsplit(rownames(x), "-")
      return(c(x[1,2:3], x[2,1], as.numeric(y[[1]][1]), as.numeric(y[[1]][2])))
    } else {
      return(c(x[1,2:3], x[2,1]))
    }
  }
  else{
    {
      if(shared.alleles){return(c(NA,NA,NA,NA,NA))} else return(c(NA,NA,NA))
    }
  }
}


#' Aggregate matrix cells (lower the resolution by a factor)
#'
#' @param void internal function to be documented
#' @keywords internal
#' @export
aggregate_matrix <- function(M, fact){
  id <- seq(1,ncol(M), by = fact)
  id <- cbind(id, c(id[-1]-1, ncol(M)))
  R <- matrix(NA, nrow(id), nrow(id))
  for(i in 1:(nrow(id)-1)){
    for(j in (i+1):nrow(id)){
      R[j,i] <-  R[i,j] <- mean(M[id[i,1]:id[i,2], id[j,1]:id[j,2]], na.rm = TRUE)
    }
  }
  R
}
