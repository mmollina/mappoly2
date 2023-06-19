#' Drop marker(s) from a sequence
#'
#' @param void internal function
#' @keywords internal
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
drop_marker <- function(input.seq,
                        mrk,
                        reestimate.map = TRUE){
  if(!mappoly2:::is.mapped.sequence(input.seq)){
    message("you can only drop markers in mapped sequences\nreturning original sequence.")
    return(input.seq)
  }
  if(is.character(mrk))
    id <- match(mrk, rownames(input.seq$phases[[1]]$p1))
  if(any(is.na(id))){
    message("'mrk' is not in the original sequence. \nreturning original sequence.")
    return(input.seq)
  }
  assert_that(is.logical(reestimate.map))
  mrk.names <- rownames(input.seq$phases[[1]]$p1)
  d <- cumsum(imf_h(c(0,input.seq$phases[[1]]$rf)))
  names(d) <- mrk.names
  if(any(!id%in%seq_along(mrk.names))){
    message("marker number(s) is not in within the range of the original sequence. \nreturning original sequence.")
    return(input.seq)
  }
  input.seq$phases[[1]]$p1 <- input.seq$phases[[1]]$p1[-id, ]
  input.seq$phases[[1]]$p2 <-input.seq$phases[[1]]$p2[-id, ]
  input.seq$phases[[1]]$loglike <- NULL
  d <- mf_h(diff(d[-id]))
  names(d) <- NULL
  input.seq$phases[[1]]$rf <- d
  input.seq$phases[[1]]$error
  input.seq$phases[[1]]$haploprob <- NULL
  if(reestimate.map)
    input.seq <- mapping(input.seq, rf = input.seq$phases[[1]]$rf)
  return(input.seq)
}

#' For a given marker, test all possible phases.
#'
#' @param void internal function
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
test_marker <- function(input.seq,
                        mrk,
                        input.mat = NULL,
                        tol = 10e-4,
                        verbose = TRUE){
  assert_that(length(mrk)==1)
  assert_that(is.haplotype.sequence(input.seq))
  assert_that(is.character(mrk))
  assert_that(mrk%in%input.seq$data$mrk.names)
  if(mrk%in%rownames(input.seq$phases[[1]]$p1)){
    message("'mrk' is alredy in the map. Retutning original sequence.")
    return(input.seq)
  }
  mrk.pos <- rownames(input.seq$phases[[1]]$p1) # positioned markers
  if(is.null(input.mat)){
    S2 <- S1 <- matrix(NA, 1, length(mrk.pos), dimnames = list(mrk, mrk.pos))
  } else {
    assert_that(matrix_contain_data_seq(input.mat,input.seq))
    S1 <- input.mat$Sh.p1[mrk, mrk.pos, drop = FALSE]
    S2 <- input.mat$Sh.p2[mrk, mrk.pos, drop = FALSE]
  }
  ## two-point phasing parent 1
  dose.vec1 <- input.seq$data$dosage.p1[mrk]
  InitPh1 <- input.seq$phases[[1]]$p1
  L1 <- mappoly2:::phasing_one(mrk, dose.vec1, S1, InitPh1, FALSE)
  ## two-point phasing parent 2
  dose.vec2 <- input.seq$data$dosage.p2[mrk]
  InitPh2 <- input.seq$phases[[1]]$p2
  L2 <- mappoly2:::phasing_one(mrk, dose.vec2, S2, InitPh2, FALSE)
  ## Selecting phase configurations
  n.conf <- sapply(L1, nrow) * sapply(L2, nrow)
  pedigree <- matrix(rep(c(1,
                           2,
                           input.seq$data$ploidy.p1,
                           input.seq$data$ploidy.p2, 1),
                         input.seq$data$n.ind),
                     nrow = input.seq$data$n.ind,
                     byrow = TRUE)
  flanking <- mappoly2:::find_flanking_markers(input.seq$mrk.names, mrk.pos, mrk)
  G <- input.seq$data$geno.dose[mrk, ,drop = TRUE]
  G[is.na(G)] <- -1
  u <- match(unlist(flanking[[mrk]]), mrk.pos)
  if(is.na(u)[1]){ # Marker inserted at the beginning of the linkage group
    homolog_prob <- as.matrix(input.seq$phases[[1]]$haploprob[,c(na.omit(u), na.omit(u)+1)+2])
    idx <- c(1,0,2)
  } else if(is.na(u)[2]){ # Marker inserted at the end of the linkage group
    homolog_prob <- as.matrix(input.seq$phases[[1]]$haploprob[,c(na.omit(u)-1, na.omit(u))+2])
    idx <- c(0,2,1)
  } else { # Marker inserted in the middle of the linkage group
    homolog_prob <- as.matrix(input.seq$phases[[1]]$haploprob[,u+2])
    idx <- c(0,1,2)
  }
  w2<-w1<-NULL
  z<-vector("list", nrow(L1[[1]]) * nrow(L2[[1]]))
  count <- 1
  for(j in 1:nrow(L1[[1]])){
    for(k in 1:nrow(L2[[1]])){
      PH <- list(L1[[1]][j,], L2[[1]][k,])
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
      w1 <- rbind(w1, L1[[1]][j,])
      w2 <- rbind(w2, L2[[1]][k,])
      count <- count + 1
    }
  }
  x <- sapply(z, function(x) x[[1]])
  x <- max(x) - x
  id <- order(x)
  phasing_results <- list(loglike = x[id],
                          rf.vec = t(sapply(z[id],
                                            function(x) x[[2]])),
                          phases = list(p1 = w1[id,,drop=FALSE],
                                        p2 = w2[id,,drop=FALSE]))
  phasing_results
}

#' Add a marker to a pre-mapped sequence
#'
#' @param void internal function
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
add_marker <- function(input.seq,
                       mrk,
                       input.mat = NULL,
                       tol = 10e-4,
                       verbose = TRUE,
                       thresh.LOD.ph.to.insert = 10,
                       thresh.rf.to.add = NULL){
  ph.res <- test_marker(input.seq,mrk,input.mat,tol,verbose)

  if(length(ph.res$loglike) != 1 & diff(ph.res$loglike[1:2]) < thresh.LOD.ph.to.insert){
    if(verbose) message("there are phase configs. with LOD < 'thresh.LOD.ph.to.insert'\nreturning original seqeunce.")
    input.seq
  }
  if(is.null(thresh.rf.to.add))
    thresh.rf.to.add <- max(input.seq$phases[[1]]$rf)
  if(thresh.rf.to.add < 0 || thresh.rf.to.add >= 0.5)
    stop("'thresh.rf.to.add' argument must be between 0 and 0.5")
  if(max(ph.res$rf.vec[1,]) > thresh.rf.to.add){
    if(verbose) message("the recombination fraction excedes 'thresh.rf.to.add'\nreturning original seqeunce.")
    input.seq
  }
  pos <- mappoly2:::find_flanking_markers(input.seq$mrk.names,
                                          rownames(input.seq$phases[[1]]$p1),
                                          mrk)
  cur.mrk <- rownames(input.seq$phases[[1]]$p1)
  if(is.na(pos[[1]]$preceding))# beginning
  {
    input.seq$phases[[1]]$p1 <- rbind(ph.res$phases$p1[1,], input.seq$phases[[1]]$p1)
    input.seq$phases[[1]]$p2 <- rbind(ph.res$phases$p2[1,], input.seq$phases[[1]]$p2)
    rownames(input.seq$phases[[1]]$p1) <- rownames(input.seq$phases[[1]]$p2) <- c(mrk, cur.mrk)
    input.seq$phases[[1]]$rf <- c(ph.res$rf.vec[1,1], input.seq$phases[[1]]$rf)
  }
  else if (is.na(pos[[1]]$succeeding)){ #end
    input.seq$phases[[1]]$p1 <- rbind(input.seq$phases[[1]]$p1, ph.res$phases$p1[1,])
    input.seq$phases[[1]]$p2 <- rbind(input.seq$phases[[1]]$p2, ph.res$phases$p2[1,])
    rownames(input.seq$phases[[1]]$p1) <- rownames(input.seq$phases[[1]]$p2) <- c(cur.mrk, mrk)
    input.seq$phases[[1]]$rf <- c(input.seq$phases[[1]]$rf, ph.res$rf.vec[1,2])
  }
  else {
    idp <- 1:match(pos[[1]]$preceding, cur.mrk)
    ids <- (match(pos[[1]]$succeeding, cur.mrk)): length(cur.mrk)
    preceding <- cur.mrk[idp]
    succeeding <- cur.mrk[ids]
    input.seq$phases[[1]]$p1 <- rbind(input.seq$phases[[1]]$p1[preceding,],
                                      ph.res$phases$p1[1,],
                                      input.seq$phases[[1]]$p1[succeeding,])
    input.seq$phases[[1]]$p2 <- rbind(input.seq$phases[[1]]$p2[preceding,],
                                      ph.res$phases$p2[1,],
                                      input.seq$phases[[1]]$p2[succeeding,])
    rownames(input.seq$phases[[1]]$p1) <- rownames(input.seq$phases[[1]]$p2) <- c(preceding, mrk, succeeding)
    input.seq$phases[[1]]$rf <- c(input.seq$phases[[1]]$rf[idp[-length(idp)]],
                                  ph.res$rf.vec[1,],
                                  input.seq$phases[[1]]$rf[ids[-length(ids)]])
  }
  return(input.seq)
}


