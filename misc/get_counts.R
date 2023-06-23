
get_counts_single_parent <- mappoly:::get_counts_one_parent


# x = c(2, 2)
# ploidy1 = 4
# ploidy2 = 2
# p.k = c(0,1)
# p.k1 = c(0,1)
# q.k = 0
# q.k1 = 0
# verbose = FALSE
# joint.prob = FALSE
# get_counts_two_parents(x, ploidy1, ploidy2, p.k,
#                       p.k1, q.k, q.k1)
get_counts_two_parents <- function(x,
                                   ploidy1,
                                   ploidy2,
                                   p.k,
                                   p.k1,
                                   q.k,
                                   q.k1,
                                   verbose = FALSE,
                                   joint.prob = FALSE) {
  #genotype progenie mrk1
  gen.prog.mk1 <- x[1]
  #genotype progenie mrk2
  gen.prog.mk2 <- x[2]
  ## Possible dosages in gametes from P
  dpk <- 0:length(p.k)
  dpk1 <- 0:length(p.k1)
  ## Possible dosages in gametes from Q
  dqk <- 0:length(q.k)
  dqk1 <- 0:length(q.k1)
  ## Combining dosages from P and Q (locus k)
  comb.all.gam.k <- expand.grid(dpk, dqk)
  ## Combining dosages from P and Q (locus k+1)
  comb.all.gam.k1 <- expand.grid(dpk1, dqk1)
  ## Combination of gametes that have x[1] doses in k
  pos.k <- comb.all.gam.k[apply(comb.all.gam.k, 1, sum)  ==  x[1], ]
  ## Combination of gametes that have x[2] doses in k+1
  pos.k1 <- comb.all.gam.k1[apply(comb.all.gam.k1, 1, sum)  ==  x[2], ]
  r <- NULL
  den <- 0
  for (i in 1:nrow(pos.k)) {
    b <- NULL
    for (j in 1:nrow(pos.k1)) {
      a1 <- get_counts_single_parent(ploidy1, p.k, p.k1, pos.k[i, 1], pos.k1[j, 1])
      a2 <- get_counts_single_parent(ploidy2, q.k, q.k1, pos.k[i, 2], pos.k1[j, 2])
      r <- rbind(r, kronecker(a1[-(2 + ploidy1/2)], a2[-(2 + ploidy2/2)]))
      b <- c(b, a1[2 + ploidy1/2] * a2[2 + ploidy2/2])
    }
    den <- den + mean(b)
  }
  r <- apply(r, 2, sum)
  y <- apply(expand.grid(0:(ploidy1/2), 0:(ploidy2/2)), 1, function(x) paste(x, collapse = ""))
  names(r) <- y
  #res <- NULL
  #for (i in sort(unique(y))) res <- c(res, sum(r[names(r) == i]))
  res <- r
  if (!joint.prob)
    res <- res/den
  #names(res) <- sort(unique(y))
  names(res) <- y
  res
}

##get_counts(4,2,c(0), c(0,1), c(0), c(0))
## ---o----o---       ---o----o---
## ---o----o---   x   ------------
## ------------
## ------------

# ploidy1 = 4
# ploidy2 = 2
# P.k = c(0,1)
# P.k1 = c(0,1)
# Q.k = 0
# Q.k1 = 0
# verbose = FALSE
# make.names = FALSE
# joint.prob = FALSE

get_counts <- function(ploidy1,
                       ploidy2,
                       P.k = NULL,
                       P.k1 = NULL,
                       Q.k = NULL,
                       Q.k1 = NULL) {

  # Adjust dP.k, dP.k1, dQ.k, and dQ.k1 to use ploidy1 and ploidy2
  if (all(is.null(P.k)))
    dP.k <- 0 else if (length(P.k) > ploidy1/2)
      dP.k <- (ploidy1/2):(ploidy1/2 + length(P.k) - ploidy1) else dP.k <- 0:length(P.k)
      if (all(is.null(P.k1)))
        dP.k1 <- 0 else if (length(P.k1) > ploidy1/2)
          dP.k1 <- (ploidy1/2):(ploidy1/2 + length(P.k1) - ploidy1) else dP.k1 <- 0:length(P.k1)

          if (all(is.null(Q.k)))
            dQ.k <- 0 else if (length(Q.k) > ploidy2/2)
              dQ.k <- (ploidy2/2):(ploidy2/2 + length(Q.k) - ploidy2) else dQ.k <- 0:length(Q.k)

              if (all(is.null(Q.k1)))
                dQ.k1 <- 0 else if (length(Q.k1) > ploidy2/2)
                  dQ.k1 <- (ploidy2/2):(ploidy2/2 + length(Q.k1) - ploidy2) else dQ.k1 <- 0:length(Q.k1)

                  # Update the get_counts_two_parents function call to use ploidy1 and ploidy2
                  counts <- NULL
                  bla <- sort(unique(kronecker(dP.k, dQ.k, "+")))
                  ble <- sort(unique(kronecker(dP.k1, dQ.k1, "+")))
                  bli <- expand.grid(ble, bla)[, 2:1]
                  blo <- bli[1:ceiling(nrow(bli)/2), ]
                  counts <- t(apply(blo, 1, get_counts_two_parents, ploidy1 = ploidy1,
                                    ploidy2 = ploidy2, p.k = P.k, p.k1 = P.k1, q.k = Q.k,
                                    q.k1 = Q.k1, joint.prob = FALSE))
                  if (nrow(bli)  ==  1) {
                    rownames(counts) <- apply(bli, 1, paste, collapse = " ")
                    return(counts)
                  }
                  if (nrow(bli)%%2  ==  1) {
                    counts <- rbind(counts, counts[(nrow(counts) - 1):1, ])
                  } else {
                    counts <- rbind(counts, counts[nrow(counts):1, ])
                  }
                  rownames(counts) <- apply(bli, 1, paste, collapse = " ")
                  return(counts)
}

# ploidy1 <- 4
# ploidy2 <- 6
# x <- c(1,1,6,5)
# get_counts_all_phases(x,ploidy1,ploidy2)
get_counts_all_phases <- function(x,
                                  ploidy1,
                                  ploidy2,
                                  verbose = FALSE,
                                  make.names = FALSE,
                                  joint.prob = FALSE) {
  pk <- x[1]
  pk1 <- x[2]
  qk <- x[3]
  qk1 <- x[4]
  if (any(is.na(c(ploidy1, ploidy2, pk, pk1, qk, qk1))))
    return(NULL)
  if (any(c(pk, pk1)  ==  0))
    sh.p <- 0 else {
      sh.p <- min(pk, pk1):0
      if (length(sh.p) > ploidy1 - max(pk, pk1))
        sh.p <- sh.p[1:(ploidy1 - max(pk, pk1) + 1)]
    }
  if (any(c(qk, qk1)  ==  0))
    sh.q <- 0 else {
      sh.q <- min(qk, qk1):0
      if (length(sh.q) > ploidy2 - max(qk, qk1))
        sh.q <- sh.q[1:(ploidy2 - max(qk, qk1) + 1)]
    }
  if (pk  ==  0)
    pk <- NULL else pk <- 0:(pk - 1)
  if (pk1  ==  0)
    pk1 <- NULL else pk1 <- 0:(pk1 - 1)
  if (qk  ==  0)
    qk <- NULL else qk <- 0:(qk - 1)
  if (qk1  ==  0)
    qk1 <- NULL else qk1 <- 0:(qk1 - 1)
  pk.ph <- NULL
  pk1.ph <- NULL
  if (length(pk) < length(pk1)) {
    for (i in 0:(length(sh.p) - 1)) {
      pk.ph <- rbind(pk.ph, pk)
      pk1.ph <- rbind(pk1.ph, pk1 + i)
    }
  } else {
    for (i in 0:(length(sh.p) - 1)) {
      if (!is.null(pk))
        pk.ph <- rbind(pk.ph, pk + i)
      pk1.ph <- rbind(pk1.ph, pk1)
    }
  }
  qk.ph <- NULL
  qk1.ph <- NULL
  if (length(qk) < length(qk1)) {
    for (i in 0:(length(sh.q) - 1)) {
      qk.ph <- rbind(qk.ph, qk)
      qk1.ph <- rbind(qk1.ph, qk1 + i)
    }
  } else {
    for (i in 0:(length(sh.q) - 1)) {
      if (!is.null(qk))
        qk.ph <- rbind(qk.ph, qk + i)
      qk1.ph <- rbind(qk1.ph, qk1)
    }
  }
  pk.num <- NULL
  if (any(is.null(pk.ph), is.null(pk1.ph)))
    pk.num <- 0 else {
      for (i in 1:nrow(pk.ph)) pk.num <- c(pk.num, sum(!is.na(match(pk.ph[i, ], pk1.ph[i, ]))))
    }
  qk.num <- NULL
  if (any(is.null(qk.ph), is.null(qk1.ph)))
    qk.num <- 0 else {
      for (i in 1:nrow(qk.ph)) qk.num <- c(qk.num, sum(!is.na(match(qk.ph[i, ], qk1.ph[i, ]))))
    }
  a.names <- expand.grid(qk.num, pk.num)
  a <- vector("list", length(pk.num) * length(qk.num))
  names(a) <- apply(a.names, 1, function(x) paste(rev(x), collapse = "-"))
  for (i in 1:length(pk.num)) {
    for (j in 1:length(qk.num)) {
      if (verbose)
        print(names(a)[(i - 1) * length(qk.num) + j])
      a[[(i - 1) * length(qk.num) + j]] <- get_counts(ploidy1,
                                                      ploidy2,
                                                      pk.ph[i, ],
                                                      pk1.ph[i, ],
                                                      qk.ph[j, ],
                                                      qk1.ph[j, ])
    }
  }
  a
}

cache_counts_twopt <- function(ploidy1, ploidy2, ncpus = 1, verbose = TRUE) {
  x1 <- as.matrix(unique(expand.grid(0:ploidy1, 0:ploidy1)))
  x2 <- as.matrix(unique(expand.grid(0:ploidy2, 0:ploidy2)))
  outm <- matrix(NA, nrow(x1)*nrow(x2), 4)
  for(j in 1:nrow(x2))
    for(i in 1:nrow(x1))
      outm[(j - 1) * nrow(x1) + i,] <- c(x1[i,], x2[j,])
  dose.names <- unique(apply(outm,1, paste, collapse = "-"))
  outm <- matrix(unlist(lapply(strsplit(dose.names, split = "-"), as.numeric)), ncol = 4, byrow = TRUE)
  dimnames(outm) <- list(paste("Conf.", 1:nrow(outm), sep = ""), c("P.k", "P.k+1", "Q.k", "Q.k+1"))
  if (verbose) {
    cat("\n   Caching the following dosage combination: \n")
    print(outm)
  }
  if (ncpus > 1) {
    if (verbose)
      cat("INFO: Using ", ncpus, " CPU's for calculation.\n")
    cl <- makeCluster(ncpus)
    on.exit(stopCluster(cl))
    y <- parApply(cl, outm, 1, get_counts_all_phases, ploidy = input.seq$ploidy, joint.prob = joint.prob)
  } else {
    if (verbose)
      cat("INFO: Going singlemode. Using one Core/CPU/PC for calculation.\n")
    y <- apply(outm, 1, get_counts_all_phases, ploidy1, ploidy2)
  }
  simple_form <- function(d,p)  abs(abs(d - (p / 2)) - (p / 2))
  x <- t(rbind(apply(outm[,1:2], 1, function(x,p) sum(as.logical(simple_form(x,p))), p = ploidy1),
               apply(outm[,3:4], 1, function(x,p) sum(as.logical(simple_form(x,p))), p = ploidy2)))
  id <- apply(x==2, 1, any)
  names(y) <- dose.names
  y[!id] <- NA
  structure(y, class = "cache.info.mappoly2")
}

print.cache.info.mappoly2 <- function(x)
  invisible(x)



a <- matrix(c(2,2,2,4,4,6,2,4,6,4,6,6), 6, 2)
full_counts <- vector("list", nrow(a))
names(full_counts) <- apply(a, 1, paste, collapse = "x")
for(i in 1:length(full_counts)){
  cat(names(full_counts)[i], "\n")
  full_counts[[i]] <- cache_counts_twopt(ploidy1 = a[i,1],
                                         ploidy2 = a[i,2],
                                         verbose = FALSE)
}
full_counts <- mappoly2:::full_counts
save(full_counts, file = "~/repos/official_repos/mappoly2/R/sysdata.rda", compress = "xz")

