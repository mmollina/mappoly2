#' Pairwise Two-Point Analysis
#'
#' Performs a pairwise two-point analysis for selected markers in a dataset.
#' The function estimates the recombination fraction for all possible linkage phase
#' configurations and their respective LOD Scores, returning the most likely one.
#'
#' @param input.data An object of class \code{mappoly2.data}
#' @param mrk.scope Specifies the range of markers for which the pairwise recombination
#'   fractions will be calculated. Acceptable values are "all", "per.chrom", and "chrom".
#'   See details for more information.
#' @param chrom Specifies the particular chromosome for which the
#'  pairwise recombination fractions will be calculated. This argument is required when
#'  the \code{mrk.scope} argument is set to \code{"chrom"}.
#'
#' @param ncpus Number of parallel processes (cores) to use (default = 1)
#' @param tol Desired accuracy. See \code{optimize()} for details
#' @param ll If \code{TRUE}, returns log-likelihood instead of LOD scores (for internal use)
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
#'   dat <- filter_data(B2721, mrk.thresh = .08, ind.thresh = 0.06)
#'   dat <- pairwise_rf(dat, mrk.scope = "chrom", chrom = "ch10")
#'
#' @export pairwise_rf
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
pairwise_rf <- function(input.data,
                        mrk.scope = c("all","per.chrom", "chrom"),
                        chrom = NULL,
                        ncpus = 1L,
                        verbose = TRUE,
                        tol = .Machine$double.eps^0.25)
{
  assert_that(mappoly2:::has.mappoly2.screened(input.data),
              msg = "The input data is not screened")
  mrk.scope <- match.arg(mrk.scope)

  if(mrk.scope == "all"){
    seq.num <- mappoly2:::get_screened_mrk_indices(input.data)
    mrk.pairs <- combn(sort(seq.num), 2)
    input.data$pairwise.rf <- pairwise_rf_full_mat(input.data = input.data,
                                                   ncpus = ncpus,
                                                   mrk.pairs = mrk.pairs,
                                                   seq.mrk.names = input.data$mrk.names[seq.num],
                                                   tol = tol,
                                                   mrk.scope = mrk.scope)
    class(input.data) <- c(class(input.data), "pairwise.rf")
    return(input.data)
  } else if (mrk.scope == "per.chrom") {
    ch <- unique(input.data$chrom[input.data$chrom != "NoChr"])
    id.num <- mappoly2:::get_mrk_indices_from_chrom(input.data, ch)
    rec.mat <- lod.mat <- lod.ph.mat <- Sh.p1 <- Sh.p2 <- matrix(NA, length(id.num), length(id.num),
                                                                 dimnames = list(input.data$mrk.names[id.num],
                                                                                 input.data$mrk.names[id.num]))
    cte <- 1
    for(i in ch){
      if(verbose) cat("  -->", i)
      seq.num <- mappoly2:::get_mrk_indices_from_chrom(input.data, i)
      mrk.pairs <- combn(sort(seq.num), 2)
      m <- mappoly2:::pairwise_rf_full_mat(input.data, ncpus, mrk.pairs,
                                           input.data$mrk.names[seq.num],
                                           tol)
      rec.mat[cte:(length(seq.num)+cte-1), cte:(length(seq.num)+cte-1)] <- m$rec.mat
      lod.mat[cte:(length(seq.num)+cte-1), cte:(length(seq.num)+cte-1)] <- m$lod.mat
      lod.ph.mat[cte:(length(seq.num)+cte-1), cte:(length(seq.num)+cte-1)] <- m$lod.ph.mat
      Sh.p1[cte:(length(seq.num)+cte-1), cte:(length(seq.num)+cte-1)] <- m$Sh.p1
      Sh.p2[cte:(length(seq.num)+cte-1), cte:(length(seq.num)+cte-1)] <- m$Sh.p2
      cte <- length(seq.num) + cte
      cat("\n")
    }
    input.data$pairwise.rf <- list(rec.mat = rec.mat,
                                   lod.mat = lod.mat,
                                   lod.ph.mat = lod.ph.mat,
                                   Sh.p1 = Sh.p1,
                                   Sh.p2 = Sh.p2,
                                   mrk.scope = mrk.scope)
    class(input.data) <- c(class(input.data), "pairwise.rf")
    return(input.data)
  } else {
    assert_that(!is.null(chrom))
    seq.num <- mappoly2:::get_mrk_indices_from_chrom(input.data, chrom)
    mrk.pairs <- combn(sort(seq.num), 2)
    input.data$pairwise.rf <- pairwise_rf_full_mat(input.data,
                                                   ncpus,
                                                   mrk.pairs,
                                                   input.data$mrk.names[seq.num],
                                                   tol,
                                                   mrk.scope)
    class(input.data) <- c(class(input.data), "pairwise.rf")
    return(input.data)
  }
}

v_2_m <- function(x, n){
  y <- base::matrix(NA, n, n)
  y[base::lower.tri(y)] <- as.numeric(x)
  y[base::upper.tri(y)] <- t(y)[base::upper.tri(y)]
  y
}

pairwise_rf_full_mat <- function(input.data,
                                 ncpus = 1L,
                                 mrk.pairs,
                                 seq.mrk.names,
                                 tol = .Machine$double.eps^0.25,
                                 mrk.scope = NULL)
{
  mrk.pairs <- mrk.pairs - 1
  count.cache <- mappoly2:::full_counts[[paste(sort(unlist(input.data[1:2])), collapse = "x")]]
  RcppParallel::setThreadOptions(numThreads = ncpus)
  count.vector = unlist(count.cache)
  count.phases = unlist(lapply(count.cache, function(x) paste0(names(x), collapse = '/')))
  count.matrix.rownames = unlist(lapply(count.cache, function(x) paste0(rownames(x[[1]]), collapse = '/')))
  count.matrix.number = unlist(lapply(count.cache, length))
  count.matrix.length = unlist(lapply(count.cache, function(x) length(c(unlist(x)) )))
  count.matrix.pos = cumsum(c(1, count.matrix.length[-length(count.matrix.length)]))
  ploidy.p1 <- input.data$ploidy.p1
  ploidy.p2 <- input.data$ploidy.p2
  if(ploidy.p1 <= ploidy.p2){
    dose.p1 <- input.data$dosage.p1
    dose.p2 <- input.data$dosage.p2
    swap.parents <- FALSE
  } else {
    ploidy.p1 <- input.data$ploidy.p2
    ploidy.p2 <- input.data$ploidy.p1
    dose.p1 <- input.data$dosage.p2
    dose.p2 <- input.data$dosage.p1
    swap.parents <- TRUE
  }
  geno <- as.matrix(input.data$geno.dose)
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
