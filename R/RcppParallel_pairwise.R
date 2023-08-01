#' @export pairwise_rf
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
pairwise_rf <- function(input.seq,
                        ncpus = 1L,
                        mrk.pairs = NULL,
                        verbose = TRUE,
                        tol = .Machine$double.eps^0.25)
{
  assert_that(is.mappoly2.sequence(input.seq))

  if (any(duplicated(input.seq$mrk.names)))
    stop("There are duplicated markers in the sequence")
  if (is.null(mrk.pairs)) {
    seq.num <- sort(mappoly2:::get_seq_indices(input.seq))
    mrk.pairs <- combn(sort(seq.num), 2) - 1
  } else {
    mrk.pairs <- mrk.pairs - 1
  }
  count.cache <- mappoly2:::full_counts[[paste(sort(unlist(input.seq$data[1:2])), collapse = "x")]]
  RcppParallel::setThreadOptions(numThreads = ncpus)
  count.vector = unlist(count.cache)
  count.phases = unlist(lapply(count.cache, function(x) paste0(names(x), collapse = '/')))
  count.matrix.rownames = unlist(lapply(count.cache, function(x) paste0(rownames(x[[1]]), collapse = '/')))
  count.matrix.number = unlist(lapply(count.cache, length))
  count.matrix.length = unlist(lapply(count.cache, function(x) length(c(unlist(x)) )))
  count.matrix.pos = cumsum(c(1, count.matrix.length[-length(count.matrix.length)]))
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
  geno <- as.matrix(input.seq$data$geno.dose[input.seq$mrk.names,])
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
                                          tol_R = tol, threads_R = ncpus)
  res[res == -1] = NA
  colnames(res) = c("Sh_P1","Sh_P2","rf","LOD_rf","LOD_ph")
  seq.mrk.names <- input.seq$mrk.names
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
  input.seq$pairwise = list(rec.mat = rf,
                            lod.mat = abs(LOD_rf),
                            lod.ph.mat = abs(LOD_ph),
                            Sh.p1 = Sh_P1,
                            Sh.p2 = Sh_P2)
  return(input.seq)
}



