#' Multi-Locus Map Estimation
#'
#' Performs multi-locus map estimation on mappoly2 data objects. This function supports
#' both sequential and parallel processing for handling large datasets efficiently.
#'
#' @param x An object representing mappoly2 data, typically containing genetic information
#'          and maps for linkage groups.
#' @param lg Optional; a vector of linkage group indices to be processed. If NULL, all
#'           linkage groups in the data object are considered.
#' @param type Character vector indicating the type of mapping to perform. Options include
#'             "mds" (multi-dimensional scaling) and "genome". Default is c("mds", "genome").
#' @param phase.conf A configuration parameter for phase determination.
#' @param rf Recombination fraction, used in the mapping calculations.
#' @param error Error tolerance in the mapping process.
#' @param verbose Logical; if TRUE, detailed progress information will be printed during
#'                processing. Default is TRUE.
#' @param tol Tolerance level for the mapping algorithm.
#' @param ret_H0 Logical; if TRUE, some hypothesis testing result is returned.
#' @param ncpus Integer; specifies the number of cores to use for parallel processing
#'
#' @return Returns an updated mappoly2 data object with mapped linkage groups.
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
#'
#' @importFrom parallel makeCluster stopCluster detectCores
#' @seealso \link[mappoly2]{mapping_one} for the function used in mapping calculations.
mapping <- function(x,
                    lg = NULL,
                    type = c("mds", "genome", "custom"),
                    phase.conf = "all",
                    rf = NULL,
                    error = 0.0,
                    ncpus = 1,
                    verbose = TRUE,
                    tol = 10e-4,
                    ret_H0 = FALSE)
{
  y <- mappoly2:::parse_lg_and_type(x,lg,type)
  g <- x$data$geno.dose
  g[is.na(g)] <- -1
  ploidy.p1 <- x$data$ploidy.p1
  ploidy.p2 <- x$data$ploidy.p2
  ind.names <- x$data$screened.data$ind.names
  cat("Multi-locus map estimation\n")

  mapData <- lapply(y$lg, function(i) {
    mrk.id <- rownames(x$maps[[i]][[type]]$phase[[1]]$p1)
    gtemp <- g[mrk.id, ind.names]
    ph <- x$maps[[i]][[y$type]]$phase
    list(gtemp = gtemp, ph = ph, ploidy.p1 = ploidy.p1, ploidy.p2 = ploidy.p2,
         phase.conf = phase.conf, rf = rf, error = error, verbose = verbose,
         tol = tol, ret_H0 = ret_H0, type = type)
  })

  mapFunc <- function(data) {
    if(data$verbose) cat("Processing linkage group\n")
    return(mappoly2:::mapping_one(data$gtemp, data$ph, data$ploidy.p1,
                                  data$ploidy.p2, data$phase.conf,
                                  "both", data$rf, data$error,
                                  data$verbose, data$tol, data$ret_H0))
  }
  if(ncpus > 1) {
    cl <- makeCluster(ncpus)
    results <- parLapply(cl, mapData, mapFunc)
    stopCluster(cl)
  } else {
    results <- lapply(mapData, mapFunc)
  }
  for(i in seq_along(y$lg))
    x$maps[[y$lg[i]]][[y$type]]$phase <- results[[i]]
  return(x)
}


mapping_one <- function(g, ph, ploidy.p1, ploidy.p2, phase.conf = "all",
                        parent.info = c("both", "p1", "p2"),
                        rf = NULL, error = 0.0, verbose = TRUE,
                        tol = 10e-4, ret_H0 = FALSE){
  if(all(phase.conf == "all"))
    phase.conf <- 1:length(ph)
  assert_that(all(phase.conf%in%1:length(ph)),
              msg = "invalid phases specified in 'phase.conf'")
  n.ind <- ncol(g)
  mrk.id <- rownames(ph[[1]]$p1)
  if(is.null(rf))
    rf <- rep(0.01, nrow(g) - 1)
  assert_that(length(rf) == nrow(g) - 1)
  if (parent.info == "both"){ ###FIXME: include detect_parent_info
    for(i in phase.conf){
      cat("   Conf.", i,":")
      pedigree <- matrix(rep(c(1,2,ploidy.p1,ploidy.p2, 1),n.ind),
                         nrow = n.ind,
                         byrow = TRUE)
      w <- mappoly2:::est_hmm_map_biallelic(PH = list(ph[[i]]$p1,
                                                      ph[[i]]$p2),
                                            G = g,
                                            pedigree = pedigree,
                                            rf = rf,
                                            err = error,
                                            verbose = verbose,
                                            detailed_verbose = FALSE,
                                            tol = tol,
                                            ret_H0 = ret_H0)
      ph[[i]]$loglike <- w[[1]]
      ph[[i]]$rf <- w[[2]]
      ph[[i]]$error <- error
    }
    cat("Done with map estimation\n")
    #return(ph)
    return(sort_phase(ph))
  }
  else if(parent.info == "p1"){
    id <- which(ploidy.p2 == dosage.p2[mrk.id])
    g[id, ] <- g[id, ] - ploidy.p2/2
    for(i in phase.conf){
      cat("   Conf.", i,":")
      w <- est_hmm_map_biallelic_single(PH = ph[[i]]$p1,
                                        G = g,
                                        rf = rf,
                                        err = error,
                                        verbose = verbose,
                                        detailed_verbose = FALSE,
                                        tol = tol,
                                        ret_H0 = ret_H0)
      ph[[i]]$loglike <- w[[1]]
      ph[[i]]$rf <- w[[2]]
      ph[[i]]$error <- error
    }
    cat("Done with map estimation\n")
    return(sort_phase(ph))
  }
  else if(parent.info == "p2") {
    id <- which(ploidy.p1 == dosage.p1[mrk.id])
    g[id, ] <- g[id, ] - ploidy.p1/2
    for(i in phase.conf){
      cat("   Conf.", i,":")
      w <- est_hmm_map_biallelic_single(PH = ph[[i]]$p2,
                                        G = g,
                                        rf = rf,
                                        err = error,
                                        verbose = verbose,
                                        detailed_verbose = FALSE,
                                        tol = tol,
                                        ret_H0 = ret_H0)
      ph[[i]]$loglike <- w[[1]]
      ph[[i]]$rf <- w[[2]]
      ph[[i]]$error <- error
    }
    cat("Done with map estimation\n")
    return(sort_phase(ph))
  }
  else {stop("it should not get here")}
}




#' Efficiently phases unprocessed markers in a given phased
#' genetic map, circumventing the need for complete
#' HMM-based map recomputation.
#'
#' @param void internal function
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
augment_phased_map <- function(x,
                               input.twopt,
                               thresh.LOD.ph = 5,
                               thresh.LOD.rf = 5,
                               thresh.rf = 0.5,
                               max.phases = 5,
                               thresh.LOD.ph.to.insert = 10,
                               thresh.rf.to.insert = NULL,
                               tol = 10e-4,
                               verbose = TRUE){
  assert_that(is.haplotype.sequence(x))
  if(all(x$mrk.names%in%rownames(x$phases[[1]]$p1))){
    message("All markers are phased. Retutning original sequence.")
    return(x)
  }
  assert_that(is.mappoly2.twopt(input.twopt))
  M <- rf_list_to_matrix(input.twopt,
                         thresh.LOD.ph = thresh.LOD.ph,
                         thresh.LOD.rf = thresh.LOD.rf,
                         thresh.rf = thresh.rf,
                         shared.alleles = TRUE)
  assert_that(matrix_contain_data_seq(M,x))
  mrk.pos <- rownames(x$phases[[1]]$p1) # positioned markers
  mrk.id <- setdiff(x$mrk.names, mrk.pos) # markers to be positioned
  ## two-point phasing parent 1
  dose.vec <- x$data$dosage.p1[mrk.id]
  InitPh1 <- x$phases[[1]]$p1
  S1 <- M$Sh.p1[mrk.id, mrk.pos]
  L1 <- mappoly2:::phasing_one(mrk.id, dose.vec, S1, InitPh1, verbose)
  ## two-point phasing parent 2
  dose.vec <- x$data$dosage.p2[mrk.id]
  InitPh2 <- x$phases[[1]]$p2
  S2 <- M$Sh.p2[mrk.id, mrk.pos]
  L2 <- mappoly2:::phasing_one(mrk.id, dose.vec, S2, InitPh2, verbose)
  ## Selecting phase configurations
  n.conf <- sapply(L1, nrow) * sapply(L2, nrow)
  if(verbose){
    cat("Distribution of phase configurations.\n")
    temp <- as.data.frame(table(n.conf))
    colnames(temp) <- c("n. phase. conf.", "frequency")
    print(temp)
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
                           x$data$ploidy.p1,
                           x$data$ploidy.p2, 1),
                         x$data$n.ind),
                     nrow = x$data$n.ind,
                     byrow = TRUE)
  flanking <- mappoly2:::find_flanking_markers(x$mrk.names, mrk.pos, mrk.id)
  phasing_results <- vector("list", length(flanking))
  names(phasing_results) <- names(flanking)
  if(verbose) pb <- txtProgressBar(min = 0, max = length(L1), style = 3)
  for(i in 1:length(L1)){
    G <- x$data$geno.dose[mrk.id[i], ,drop = TRUE]
    G[is.na(G)] <- -1
    u <- match(unlist(flanking[[mrk.id[i]]]), mrk.pos)
    if(is.na(u)[1]){ # Marker inserted at the beginning of the linkage group
      homolog_prob <- as.matrix(x$phases[[1]]$haploprob[,c(na.omit(u), na.omit(u)+1)+2])
      idx <- c(1,0,2)
    } else if(is.na(u)[2]){ # Marker inserted at the end of the linkage group
      homolog_prob <- as.matrix(x$phases[[1]]$haploprob[,c(na.omit(u)-1, na.omit(u))+2])
      idx <- c(0,2,1)
    } else { # Marker inserted in the middle of the linkage group
      homolog_prob <- as.matrix(x$phases[[1]]$haploprob[,u+2])
      idx <- c(0,1,2)
    }
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
                                                                   idx,
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
                                         phases = list(p1 = w1[id,,drop=FALSE],
                                                       p2 = w2[id,,drop=FALSE]))
    if(verbose) setTxtProgressBar(pb, i)

  }
  if(verbose) close(pb)
  selected.list<-phasing_results[sapply(phasing_results,
                                        function(x) length(x$loglike)==1 ||
                                          x$loglike[2] > thresh.LOD.ph.to.insert)]
  if(is.null(thresh.rf.to.insert))
    thresh.rf.to.insert <- max(x$phases[[1]]$rf)
  if(thresh.rf.to.insert < 0 || thresh.rf.to.insert >= 0.5)
    stop("'thresh.rf.to.insert' parameter must be between 0 and 0.5")
  selected.list <- phasing_results[sapply(phasing_results, function(x) max(x$rf.vec[1,]) <= thresh.rf.to.insert)]
  for(i in names(selected.list)){
    pos <- mappoly2:::find_flanking_markers(x$mrk.names,
                                            rownames(x$phases[[1]]$p1),
                                            i)
    if(length(unlist(pos)) == 0) next()
    cur.mrk <- rownames(x$phases[[1]]$p1)
    if(is.na(pos[[1]]$preceding))# beginning
    {
      x$phases[[1]]$p1 <- rbind(selected.list[[i]]$phases$p1[1,], x$phases[[1]]$p1)
      x$phases[[1]]$p2 <- rbind(selected.list[[i]]$phases$p2[1,], x$phases[[1]]$p2)
      rownames(x$phases[[1]]$p1) <- rownames(x$phases[[1]]$p2) <- c(i, cur.mrk)
    }
    else if (is.na(pos[[1]]$succeeding)){ #end
      x$phases[[1]]$p1 <- rbind(x$phases[[1]]$p1, selected.list[[i]]$phases$p1[1,])
      x$phases[[1]]$p2 <- rbind(x$phases[[1]]$p2, selected.list[[i]]$phases$p2[1,])
      rownames(x$phases[[1]]$p1) <- rownames(x$phases[[1]]$p2) <- c(cur.mrk, i)
    }
    else {
      preceding <- cur.mrk[1:match(pos[[1]]$preceding, cur.mrk)]
      succeeding <- cur.mrk[(match(pos[[1]]$succeeding, cur.mrk)): length(cur.mrk)]
      x$phases[[1]]$p1 <- rbind(x$phases[[1]]$p1[preceding,],
                                selected.list[[i]]$phases$p1[1,],
                                x$phases[[1]]$p1[succeeding,])
      x$phases[[1]]$p2 <- rbind(x$phases[[1]]$p2[preceding,],
                                selected.list[[i]]$phases$p2[1,],
                                x$phases[[1]]$p2[succeeding,])
      rownames(x$phases[[1]]$p1) <- rownames(x$phases[[1]]$p2) <- c(preceding, i, succeeding)
    }
  }
  return(x)
}

#'This function merges two genetic maps built from markers that are exclusively
#'informative in isolated parents, facilitating unified analysis and visualization
#'of distinct genetic data.
#'
#' @param void internal function
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
merge_single_parent_maps <- function(x.all,
                                     x.p1,
                                     x.p2,
                                     input.twopt){
  temp.mrk.id <- character(length(x.all$mrk.names))
  temp.mrk.id[match(rownames(x.p1$phases[[1]]$p1), x.all$mrk.names)] <- rownames(x.p1$phases[[1]]$p1)
  temp.mrk.id[match(rownames(x.p2$phases[[1]]$p2), x.all$mrk.names)] <- rownames(x.p2$phases[[1]]$p2)
  temp.mrk.id <- temp.mrk.id[temp.mrk.id != ""]
  ph.p2 <- ph.p1 <- matrix(0, length(temp.mrk.id), x.p1$data$ploidy.p1, dimnames = list(temp.mrk.id, NULL))
  ph.p1[rownames(x.p1$phases[[1]]$p1), ] <- x.p1$phases[[1]]$p1
  ph.p2[rownames(x.p1$phases[[1]]$p2), ] <- x.p1$phases[[1]]$p2
  ph.p1[rownames(x.p2$phases[[1]]$p1), ] <- x.p2$phases[[1]]$p1
  ph.p2[rownames(x.p2$phases[[1]]$p2), ] <- x.p2$phases[[1]]$p2
  x.all$phases <- list(list(p1 = ph.p1,
                            p2 = ph.p2,
                            loglike = NULL,
                            rf = NULL,
                            error = NULL,
                            haploprob = NULL))
  return(x.all)
}

