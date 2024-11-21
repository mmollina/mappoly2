#' Pairwise Two-Point Analysis
#'
#' Performs a pairwise two-point analysis for selected markers in a dataset.
#' The function estimates the recombination fraction for all possible linkage phase
#' configurations and their respective LOD Scores, returning the most likely one.
#'
#' @param x An object of class \code{mappoly2.data}
#' @param mrk.scope Specifies the range of markers for which the pairwise recombination
#'   fractions will be calculated. Acceptable values are "all", "per.chrom", and "chrom".
#'   See details for more information.
#' @param chrom Specifies the particular chromosome for which the
#'  pairwise recombination fractions will be calculated. This argument is required when
#'  the \code{mrk.scope} argument is set to \code{"chrom"}.
#' @param ncpus Number of parallel processes (cores) to use (default = 1)
#' @param verbose Logical; if TRUE, progress messages will be printed.
#' @param tol Desired accuracy. See \code{optimize()} for details
#'
#' @details
#' The \code{mrk.scope} argument allows for the customization of the analysis scope. Options are:
#'   - \code{"all"}: Analyze all markers in the dataset.
#'   - \code{"per.chrom"}: Analyze markers within each chromosome.
#'   - \code{"chrom"}: Analyze markers in a specific chromosome. Requires specifying the chromosome.
#'
#' Additional arguments (\code{chrom} or \code{sequence}) must be specified when using
#' \code{"chrom"} or \code{"seq"} options, respectively.
#'
#' @return An updated object of class \code{mappoly2.data} with a slot called 'pairwise.rf' containing:
#'
#'  - \code{"rec.mat"} the recombination fraction matrix
#'  - \code{"lod.mat"} the LOD Scrore associated to the recombination fraction matrix
#'  - \code{"lod.ph.mat"} the LOD Scrore associated to the second most likely linkage phase configuration
#'  - \code{"Sh.p1" and "Sh.p2"} the number of homologs that share the alternate alleles for the estimated linkage phase configurationns for parents 1 and 2, respectively.
#' @examples
#'
#'   dat <- subset(B2721, perc = .1)
#'   dat <- filter_data(dat, mrk.thresh = .08, ind.thresh = 0.06)
#'   dat <- pairwise_rf(dat, mrk.scope = "chrom", chrom = "ch10")
#'
#' @export pairwise_rf
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom graphics lines
#' @importFrom utils combn tail
pairwise_rf <- function(x,
                        mrk.scope = c("all","per.chrom", "chrom"),
                        chrom = NULL,
                        ncpus = 1L,
                        verbose = TRUE,
                        tol = .Machine$double.eps^0.25)
{
  assert_that(has.mappoly2.screened(x),
              msg = "The input data is not screened")
  mrk.scope <- match.arg(mrk.scope)
  if(!is.null(chrom))
    mrk.scope = "chrom"
  if(mrk.scope == "all"){
    seq.num <- get_screened_mrk_indices(x)
    mrk.pairs <- combn(sort(seq.num), 2)
    x$pairwise.rf <- pairwise_rf_full_mat(x = x,
                                          ncpus = ncpus,
                                          mrk.pairs = mrk.pairs,
                                          seq.mrk.names = x$mrk.names[seq.num],
                                          tol = tol,
                                          mrk.scope = mrk.scope)
    class(x) <- unique(c(class(x), "pairwise.rf"))
    return(x)
  } else if (mrk.scope == "per.chrom") {
    ch <- unique(x$chrom[x$chrom != "NoChr"])
    ch <- names(which(table(x$chrom) > 1))
    ch <- ch[order(mappoly2:::embedded_to_numeric(ch))]
    id.num <- lapply(ch, function(y) mappoly2:::get_mrk_indices_from_chrom(x, y))
    names(id.num) <- ch
    id.num <- id.num[sapply(id.num, length) > 1]
    ch <- names(id.num)
    v <- unlist(id.num)
    rec.mat <- lod.mat <- lod.ph.mat <- Sh.p1 <- Sh.p2 <- matrix(NA, length(v), length(v),
                                                                 dimnames = list(x$mrk.names[v],
                                                                                 x$mrk.names[v]))
    cte <- 1
    for(i in ch){
      if(verbose) cat("  -->", i)
      seq.num <- id.num[[i]]
      mrk.pairs <- combn(sort(seq.num), 2)
      m <- pairwise_rf_full_mat(x, ncpus, mrk.pairs,
                                x$mrk.names[seq.num],
                                tol)
      rec.mat[cte:(length(seq.num)+cte-1), cte:(length(seq.num)+cte-1)] <- m$rec.mat
      lod.mat[cte:(length(seq.num)+cte-1), cte:(length(seq.num)+cte-1)] <- m$lod.mat
      lod.ph.mat[cte:(length(seq.num)+cte-1), cte:(length(seq.num)+cte-1)] <- m$lod.ph.mat
      Sh.p1[cte:(length(seq.num)+cte-1), cte:(length(seq.num)+cte-1)] <- m$Sh.p1
      Sh.p2[cte:(length(seq.num)+cte-1), cte:(length(seq.num)+cte-1)] <- m$Sh.p2
      cte <- length(seq.num) + cte
      cat("\n")
    }
    x$pairwise.rf <- list(rec.mat = rec.mat,
                          lod.mat = lod.mat,
                          lod.ph.mat = lod.ph.mat,
                          Sh.p1 = Sh.p1,
                          Sh.p2 = Sh.p2,
                          mrk.scope = mrk.scope)
    class(x) <- unique(c(class(x), "pairwise.rf"))
    return(x)
  } else {
    assert_that(!is.null(chrom))
    seq.num <- get_mrk_indices_from_chrom(x, chrom)
    mrk.pairs <- combn(sort(seq.num), 2)
    x$pairwise.rf <- pairwise_rf_full_mat(x,
                                          ncpus,
                                          mrk.pairs,
                                          x$mrk.names[seq.num],
                                          tol,
                                          mrk.scope)
    class(x) <- unique(c(class(x), "pairwise.rf"))
    return(x)
  }
}

v_2_m <- function(x, n){
  y <- base::matrix(NA, n, n)
  y[base::lower.tri(y)] <- as.numeric(x)
  y[base::upper.tri(y)] <- t(y)[base::upper.tri(y)]
  y
}

pairwise_rf_full_mat <- function(x,
                                 ncpus = 1L,
                                 mrk.pairs,
                                 seq.mrk.names,
                                 tol = .Machine$double.eps^0.25,
                                 mrk.scope = NULL)
{
  mrk.pairs <- mrk.pairs - 1
  count.cache <- full_counts[[paste(sort(unlist(x[1:2])), collapse = "x")]]
  RcppParallel::setThreadOptions(numThreads = ncpus)
  count.vector = unlist(count.cache)
  count.phases = unlist(lapply(count.cache, function(x) paste0(names(x), collapse = '/')))
  count.matrix.rownames = unlist(lapply(count.cache, function(x) paste0(rownames(x[[1]]), collapse = '/')))
  count.matrix.number = unlist(lapply(count.cache, length))
  count.matrix.length = unlist(lapply(count.cache, function(x) length(c(unlist(x)) )))
  count.matrix.pos = cumsum(c(1, count.matrix.length[-length(count.matrix.length)]))
  ploidy.p1 <- x$ploidy.p1
  ploidy.p2 <- x$ploidy.p2
  if(ploidy.p1 <= ploidy.p2){
    dose.p1 <- x$dosage.p1
    dose.p2 <- x$dosage.p2
    swap.parents <- FALSE
  } else {
    ploidy.p1 <- x$ploidy.p2
    ploidy.p2 <- x$ploidy.p1
    dose.p1 <- x$dosage.p2
    dose.p2 <- x$dosage.p1
    swap.parents <- TRUE
  }
  geno <- as.matrix(x$geno.dose)
  geno[is.na(geno)] <- 1 + (ploidy.p1 + ploidy.p2)/2
  res <- pairwise_rf_estimation_disc_rcpp(mrk_pairs_R = as.matrix(mrk.pairs),
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
                                          tol_R = tol)
  res[res == -1] = NA
  colnames(res) = c("Sh_P1","Sh_P2","rf","LOD_rf","LOD_ph")
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
  pairwise.rf = list(rec.mat = rf,
                     lod.mat = abs(LOD_rf),
                     lod.ph.mat = abs(LOD_ph),
                     Sh.p1 = Sh_P1,
                     Sh.p2 = Sh_P2,
                     mrk.scope = mrk.scope)
  return(pairwise.rf)
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
