z <- x$working.sequences
n.seq <- length(z)
R1 <- NULL
for(i in 1:n.seq){
  n.mrk <- length(z[[i]]$mrk.names)
  n.ind <- length(z[[i]]$ind.names)
  n.mrk.mds <- nrow(z[[i]]$order$mds$info$locimap)
  n.conf.2pt.phase.mds <- length(z[[i]]$order$mds$phase.twopt)
  n.mrk.2pt.phase.mds <- nrow(z[[i]]$order$mds$phase.twopt[[1]]$p1)
  n.mrk.hmm.phase.mds <- nrow(z[[i]]$order$mds$hmm.map$p1)
  chrom.mds <- unique(x$data$chrom[rownames(z[[i]]$order$mds$hmm.map$p1)])
  map.length.mds <- round(sum(imf_h(z[[i]]$order$mds$hmm.map$rf)),1)
  max.gap.mds <- round(max(imf_h(z[[i]]$order$mds$hmm.map$rf)),1)
  a <- get_dosage_type(x, i, "mds")
  simplex.p1.mds <- length(a$simplex.p1)
  simplex.p2.mds <- length(a$simplex.p2)
  double.simplex.mds <- length(a$double.simplex)
  multiplex.mds <- length(a$multiplex)
  mrk.per.cm.mds <- round(n.mrk.hmm.phase.mds/map.length.mds,1)
  r <- c(i, chrom.mds, n.mrk, n.mrk.mds, n.conf.2pt.phase.mds, n.mrk.2pt.phase.mds,
         n.mrk.hmm.phase.mds,  map.length.mds, mrk.per.cm.mds, max.gap.mds,
         simplex.p1.mds, simplex.p2.mds, double.simplex.mds,
         multiplex.mds)
  names(r) <- c("LG", "Chrom", "Alloc Mrk", "Ordered Mrk", "N phases 2pt", "N mkr phased 2pt",
                "N mrk phased HMM", "Map Length (cM)", "Mrk/cM", "Max Gap", "Simplex P1", "Simplex P2",
                "Double Simplex", "Multiplex")
  R1 <- rbind(R1,r)
}
as.data.frame(t(R))


R2 <- NULL
for(i in 1:n.seq){
  n.mrk <- length(z[[i]]$mrk.names)
  n.ind <- length(z[[i]]$ind.names)
  n.mrk.genome <- nrow(z[[i]]$order$genome$info)
  n.conf.2pt.phase.genome <- length(z[[i]]$order$genome$phase.twopt)
  n.mrk.2pt.phase.genome <- nrow(z[[i]]$order$genome$phase.twopt[[1]]$p1)
  n.mrk.hmm.phase.genome <- nrow(z[[i]]$order$genome$hmm.map$p1)
  chrom.genome <- unique(x$data$chrom[rownames(z[[i]]$order$genome$hmm.map$p1)])
  map.length.genome <- round(sum(imf_h(z[[i]]$order$genome$hmm.map$rf)),1)
  max.gap.genome <- round(max(imf_h(z[[i]]$order$genome$hmm.map$rf)),1)
  a <- get_dosage_type(x, i, "genome")
  simplex.p1.genome <- length(a$simplex.p1)
  simplex.p2.genome <- length(a$simplex.p2)
  double.simplex.genome <- length(a$double.simplex)
  multiplex.genome <- length(a$multiplex)
  mrk.per.cm.genome <- round(n.mrk.hmm.phase.genome/map.length.genome,1)
  r <- c(i, chrom.genome, n.mrk, n.mrk.genome, n.conf.2pt.phase.genome, n.mrk.2pt.phase.genome,
         n.mrk.hmm.phase.genome,  map.length.genome, mrk.per.cm.genome, max.gap.genome,
         simplex.p1.genome, simplex.p2.genome, double.simplex.genome,
         multiplex.genome)
  names(r) <- c("LG", "Chrom", "Alloc Mrk", "Ordered Mrk", "N phases 2pt", "N mkr phased 2pt",
                "N mrk phased HMM", "Map Length (cM)", "Mrk/cM", "Max Gap", "Simplex P1", "Simplex P2",
                "Double Simplex", "Multiplex")
  R2 <- rbind(R2,r)
}

cbind(as.data.frame(t(R1)), "|", as.data.frame(t(R2)))

