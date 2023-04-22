sim_homologous <- function(ploidy1, ploidy2, n.mrk, min.d = 0,
                           max.d1 = ploidy1, max.d2 = ploidy2,
                           max.ph, restriction = TRUE,
                           seed = NULL){
  #prob.dose <- NULL
  if(!is.null(seed)) set.seed(seed)

  hom.allele.q <- hom.allele.p <- vector("list", n.mrk)
  count <- 1
  while(count <= n.mrk)
  {
    hom.p.temp <- hom.q.temp <- 0
    p.temp <- sample(min.d:max.d1,1)
    if(all(p.temp != 0))
      hom.p.temp <- sample(1:ploidy1, p.temp)
    q.temp <- sample(min.d:max.d2,1)
    if(all(q.temp != 0))
      hom.q.temp <- sample(1:ploidy2, q.temp)
    p <- sum(as.logical(hom.p.temp))
    q <- sum(as.logical(hom.q.temp))
    if(restriction && count > 1)
    {
      if(!any((p+q) == 0,
              (p+q) == ploidy1 + ploidy2,
              sum(as.logical(hom.allele.p[[count-1]])) - q  ==  0,
              sum(as.logical(hom.allele.q[[count-1]])) - p  ==  0,
              (p == ploidy1 & q == 0),
              (p == 0 & q == ploidy2)))
      {
        hom.allele.p[[count]] <- hom.p.temp
        hom.allele.q[[count]] <- hom.q.temp
        count <- count+1
      }
    }
    else
    {
      if(!any((p+q) == 0,
              (p+q) == ploidy1 + ploidy2,
              (p == ploidy1 & q == 0),
              (p == 0 & q == ploidy2)))
      {
        hom.allele.p[[count]] <- hom.p.temp
        hom.allele.q[[count]] <- hom.q.temp
        count <- count+1
      }
    }
    p <- unlist(lapply(hom.allele.p, function(x) sum(as.logical(x))))
    q <- unlist(lapply(hom.allele.q, function(x) sum(as.logical(x))))
  }
  return(list(hom.allele.p = hom.allele.p,
              hom.allele.q = hom.allele.q,
              p = p, q = q))
}
sim_cross_one_informative_parent <- function(ploidy,
                                             n.mrk,
                                             rf.vec,
                                             hom.allele,
                                             n.ind,
                                             seed = NULL,
                                             prob = NULL){
  if(!is.null(seed)) set.seed(seed)
  if(length(rf.vec) == 1) rf.vec <- rep(rf.vec, n.mrk-1)
  res <- matrix(NA,n.ind,n.mrk)
  rf.res <- numeric(n.mrk-1)

  ## Listing all possible bivalent configurations
  a <- mappoly:::perm_tot(1:ploidy)
  bv.conf <- vector("list", nrow(a))
  for(i in 1:nrow(a))
  {
    temp <- apply(matrix(a[i,], 2, ploidy/2), 2, sort)
    bv.conf[[i]] <- temp[,order(temp[1,]), drop = FALSE]
  }
  bv.conf <- unique(bv.conf)
  names(bv.conf) <- sapply(bv.conf, function(x) paste(apply(x, 2,
                                                            function(x)
                                                              paste0("[", paste0(x, collapse = ""), "]", collapse = "")),
                                                      collapse = ""))
  if(is.null(prob))
    prob <- rep(1/length(bv.conf), length(bv.conf))

  for(k in 1:n.ind){              #for each individual
    gen.1 <- matrix(1:ploidy,ploidy,n.mrk)  #simulates the chromosomes multiallelic markers in 'n.mrk' positions
    id <- sample(x = 1:length(bv.conf), size = 1, prob = prob) #sampling one bivalent configuration based on given probabilities
    choosed_biv <- bv.conf[[id]]
    choosed_biv <- choosed_biv[,sample(1:(ploidy/2)), drop = FALSE]
    for(i in 1:ncol(choosed_biv))
    {
      choosed_biv[,i] <- sample(choosed_biv[,i, drop = FALSE])
    }
    pole.1 <- choosed_biv[1,, drop = FALSE]
    pole.2 <- choosed_biv[2,, drop = FALSE]
    set.2 <- gen.1[pole.1,, drop = FALSE]      #allocating the chromosomes on the variables set.1 and set.2, thus (set.1[i], set.2[i]) represents a bivalent
    set.1 <- gen.1[pole.2,, drop = FALSE]
    for(i in 1:(ploidy/2)){         #for each one of the ploidy/2 chromosome pair (bivalents)
      a <- set.1[i,]
      b <- set.2[i,]
      for(j in 1:(n.mrk-1)){             #for each adjacent interval between.mrkkers
        if(runif(1)  < rf.vec[j]){       #if a random number drawn from the interval [0,1] (according a uniform distribution)
          #is less than the recombination fraction for that interval
          which.swap <- c((j+1):n.mrk)     #the alleles for that interval and bivalent are swapped
          temp <- a[which.swap]
          a[which.swap] <- b[which.swap]
          b[which.swap] <- temp
        }
      }             #this completes the whole bivalent
      set.1[i,] <- a  #attributing the resulting vector to the initial variables
      set.2[i,] <- b
    }               #for all bivalents
    if(sample(0:1,1)) gam <- set.1 #sample one of the meiotic products
    else gam <- set.2
    for(i in 1:(ploidy/2)){ #counting the recombinant chromosomes in their multiallelic form
      for(j in 2:ncol(gam)){
        if(!gam[i,j] == gam[i,j-1])
          rf.res[j-1] <- rf.res[j-1]+1
      }
    }
    for(i in 1:n.mrk)
      gam[,i] <- as.numeric(!is.na(match(gam[,i], hom.allele[[i]])))
    res[k,] <- apply(gam, 2, sum)
  }
  rf.calc <- rf.res/(n.ind*ploidy/2)  #computing the recombination fraction
  dimnames(res) <- list(paste("Ind",1:n.ind, sep = "_"), paste("M",1:n.mrk, sep = "_"))
  list(data.sim.one.parent = res, cross.count.one.parent = rf.res, rf.calc.one.parent = rf.calc)
}
sim_cross_two_informative_parents <- function(ploidy1,
                                              ploidy2,
                                              n.mrk,
                                              rf.vec,
                                              n.ind,
                                              hom.allele.p,
                                              hom.allele.q,
                                              prob.P = NULL,
                                              prob.Q = NULL,
                                              seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  dose.p <- unlist(lapply(hom.allele.p, function(x) sum(as.logical(x))))
  dose.q <- unlist(lapply(hom.allele.q, function(x) sum(as.logical(x))))
  if(any(apply(rbind(dose.p,dose.q),2,sum) == 0)) stop("Found zero doses in both parents at the same marker")
  if(!is.null(seed)) set.seed(seed)
  ##Parent 1
  data.P <- sim_cross_one_informative_parent(ploidy = ploidy1, n.mrk = n.mrk, rf.vec = rf.vec,
                                             hom.allele = hom.allele.p, n.ind = n.ind)
  ##Parent 2
  data.Q <- sim_cross_one_informative_parent(ploidy = ploidy2, n.mrk = n.mrk, rf.vec = rf.vec,
                                             hom.allele = hom.allele.q, n.ind = n.ind)
  rf.calc <- (data.P$cross.count.one.parent + data.Q$cross.count.one.parent)/(n.ind*(ploidy1/2 + ploidy2/2))
  data.sim.two.parents <- data.P$data.sim.one.parent + data.Q$data.sim.one.parent
  list(data.sim.two.parents = data.sim.two.parents, rf.calc = rf.calc)
}

test_simulate <- function(ploidy.p1,
                          ploidy.p2,
                          fpath = "fake_triploid.csv",
                          n.mrk = 100,
                          n.ind = 200,
                          map.length = 100){
  h <- sim_homologous(ploidy.p1, ploidy.p2, n.mrk)
  rf.vec <- mappoly::mf_h(diff(seq(0, map.length, length.out = n.mrk)))
  dat <- t(sim_cross_two_informative_parents(ploidy1 = ploidy.p1,
                                           ploidy2 = ploidy.p2,
                                           rf.vec = rf.vec,
                                           hom.allele.p = h$hom.allele.p,
                                           hom.allele.q = h$hom.allele.q,
                                           n.mrk = n.mrk, n.ind = n.ind)$data.sim.two.parents)
  x<-data.frame(snp_id = rownames(dat),
                P1 = h$p,
                P2 = h$q,
                chrom = 1,
                genome_pos = 1:n.mrk,
                ref = B2721$ref[1:n.mrk],
                alt = B2721$alt[1:n.mrk],
                dat)
  write.csv(x, fpath, row.names = FALSE)
  invisible(list(p1 = mappoly::ph_list_to_matrix(h$hom.allele.p, ploidy.p1),
              p2 = mappoly::ph_list_to_matrix(h$hom.allele.q, ploidy.p2)))
}











