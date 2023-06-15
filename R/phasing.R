#' Phasing based on pairwise recombination fraction estimation
#'
#' @param void internal function
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
pairwise_phasing <- function(input.seq,
                             input.twopt,
                             thresh.LOD.ph = 3,
                             thresh.LOD.rf = 3,
                             thresh.rf = 0.5,
                             max.conf.btnk.p1 = 1,
                             max.conf.btnk.p2 = max.conf.btnk.p1,
                             verbose = TRUE){
  m <- rf_list_to_matrix(tpt,
                         thresh.LOD.ph = thresh.LOD.ph,
                         thresh.LOD.rf = thresh.LOD.rf,
                         thresh.rf = thresh.rf,
                         shared.alleles = TRUE,
                         verbose = FALSE)
  mrk.id <- input.seq$mrk.names
  {
    if(verbose)
      cat("Phasing parent", input.seq$data$name.p1, "\n")
    Ph.p1 <- mappoly2:::twopt_phasing_cpp(mrk_id = mrk.id,
                                          ploidy = input.seq$data$ploidy.p1,
                                          dose_vec = input.seq$data$dosage.p1,
                                          S = m$Sh.p1[mrk.id, mrk.id],
                                          max_conf_number = max.conf.btnk.p1,
                                          verbose = verbose)
    for(i in 1:length(Ph.p1$phase_configs))
      rownames(Ph.p1$phase_configs[[i]]) <- Ph.p1$marker_names
    if(verbose)
      cat("Phasing parent", input.seq$data$name.p2, "\n")
    Ph.p2 <- mappoly2:::twopt_phasing_cpp(mrk_id = mrk.id,
                                          ploidy = input.seq$data$ploidy.p2,
                                          dose_vec = input.seq$data$dosage.p2,
                                          S = m$Sh.p2[mrk.id, mrk.id],
                                          max_conf_number = max.conf.btnk.p2,
                                          verbose = verbose)
    for(i in 1:length(Ph.p2$phase_configs))
      rownames(Ph.p2$phase_configs[[i]]) <- Ph.p2$marker_names
  }

  mrks <- intersect(Ph.p1$marker_names, Ph.p2$marker_names)
  n1 <- length(input.seq$mrk.names)
  n2 <- length(mrks)
  if(verbose){
    cat(n2, " phased markers out of ", n1, ": (",  round(100*n2/n1,1), "%)",sep = "")
  }
  cte <- 1
  Ph <- vector("list", length(Ph.p1$phase_configs) * length(Ph.p2$phase_configs))
  for(i in 1:length(Ph.p1$phase_configs)){
    for(j in 1:length(Ph.p2$phase_configs)){
      Ph[[cte]] <- list(p1 = Ph.p1$phase_configs[[i]][mrks,],
                        p2 = Ph.p2$phase_configs[[j]][mrks,],
                        loglike = NULL,
                        rf = NULL,
                        error = NULL,
                        haploprob = NULL)
      cte <- cte + 1
    }
  }
  input.seq$phases = unique(Ph)
  return(input.seq)
}


#' Phasing remaining markers based on pairwise recombination fraction estimation and multilocus estimation
#'
#' @param void internal function
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
phase_remaining <- function(input.seq,
                            input.twopt,
                            thresh.LOD.ph = 3,
                            thresh.LOD.rf = 3,
                            thresh.rf = 0.5,
                            max.phases = 3,
                            tol = 10e-4,
                            verbose = TRUE){
  assert_that(is.haplotype.sequence(input.seq))
  assert_that(is.mappoly2.twopt(input.twopt))
  M <- rf_list_to_matrix(input.twopt,
                         thresh.LOD.ph = thresh.LOD.ph,
                         thresh.LOD.rf = thresh.LOD.rf,
                         thresh.rf = thresh.rf,
                         shared.alleles = TRUE)
  assert_that(matrix_contain_data_seq(M,s1))
  mrk.pos <- rownames(input.seq$phases[[1]]$p1) # positioned markers
  mrk.id <- setdiff(input.seq$mrk.names, mrk.pos) # markers to be positioned
  ## two-point phasing parent 1
  dose.vec <- input.seq$data$dosage.p1[mrk.id]
  InitPh1 <- input.seq$phases[[1]]$p1
  S1 <- M$Sh.p1[mrk.id, mrk.pos]
  L1 <- mappoly2:::phasing_one(mrk.id, dose.vec, S1, InitPh1, verbose)
  ## two-point phasing parent 2
  dose.vec <- input.seq$data$dosage.p2[mrk.id]
  InitPh2 <- input.seq$phases[[1]]$p2
  S2 <- M$Sh.p2[mrk.id, mrk.pos]
  L2 <- mappoly2:::phasing_one(mrk.id, dose.vec, S2, InitPh2, verbose)
  ## Selecting phase configurations
  n.conf <- sapply(L1, nrow) + sapply(L2, nrow)
  if(verbose){
    cat("Distribution of phase configurations.\n")
    txtplot::txtboxplot(n.conf)
  }
  mrk.sel <- which(n.conf <= max.phases)
  if(length(mrk.sel) == 0)
    stop("No markers were selected for 'max.phases' = ", max.phases,
         "\n'max.phases' should be at least ", min(n.conf))
  L1 <- L1[n.conf <= max.phases]
  L2 <- L2[n.conf <= max.phases]
  mrk.id <- mrk.id[n.conf <= max.phases]
  pedigree <- matrix(rep(c(1,
                           2,
                           input.seq$data$ploidy.p1,
                           input.seq$data$ploidy.p2, 1),
                         input.seq$data$n.ind),
                     nrow = input.seq$data$n.ind,
                     byrow = TRUE)
  flanking <- find_flanking_markers(input.seq$mrk.names, mrk.pos, mrk.id)
  phasing_results <- vector("list", length(flanking))
  names(phasing_results) <- names(flanking)
  if(verbose) pb <- txtProgressBar(min = 0, max = length(L1), style = 3)

  for(i in 1:length(L1)){
    G <- input.seq$data$geno.dose[mrk.id[i], ,drop = TRUE]
    G[is.na(G)] <- -1
    u <- match(unlist(flanking[[mrk.id[i]]]), mrk.pos)
    homolog_prob <- as.matrix(input.seq$phases[[1]]$haploprob[,u+2])
    w2<-w1<-NULL
    z<-vector("list", nrow(L1[[i]]) * nrow(L2[[i]]))
    count <- 1
    for(j in 1:nrow(L1[[i]])){
      for(k in 1:nrow(L2[[i]])){
        PH <- list(L1[[i]][j,], L2[[i]][k,])
        z[[count]]<-mappoly2:::est_hmm_map_biallelic_insert_marker(PH,
                                                                   G,
                                                                   pedigree,
                                                                   homolog_prob,
                                                                   rf = c(0.01,0.01),
                                                                   verbose = FALSE,
                                                                   detailed_verbose = FALSE,
                                                                   tol = tol,
                                                                   ret_H0 = FALSE)
        w1 <- rbind(w1, L1[[i]][j,])
        w2 <- rbind(w2, L2[[i]][k,])
        count <- count + 1
      }
    }
    x <- sapply(z, function(x) x[[1]])
    x <- max(x) - x
    id <- order(x)
    phasing_results[[mrk.id[i]]] <- list(loglike = x[id],
                                         rf.vec = t(sapply(z[id],
                                                           function(x) x[[2]])),
                                         phases = list(p1 = w1[id,],
                                                       p2 = w2[id,]))
    if(verbose) setTxtProgressBar(pb, i)
  }
  if(verbose) close(pb)
  return(phasing_results)
}


