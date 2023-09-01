###### OLD ####
rm(list = ls())
require(mappoly2)
source("misc/simulation.R")
ploidy.p1 = 4
ploidy.p2 = 4
n.mrk <- 100
map.length = 30
ph<-test_simulate(ploidy.p1 = ploidy.p1,
                  ploidy.p2 = ploidy.p2,
                  fpath = "~/repos/official_repos/misc/fake_triploid.csv",
                  n.mrk = n.mrk,
                  n.ind = 300,
                  map.length = map.length,
                  miss.perc = 0,
                  n.chrom = 1,
                  random = FALSE,
                  seed = 2986)
dat <- read_geno_csv(file.in = "misc/fake_triploid.csv",
                     ploidy.p1 = ploidy.p1, ploidy.p2 = ploidy.p2)
s <- make_sequence(dat, "all", info.parent = "p1")
tpt <- est_pairwise_rf(s)
s.phased <- pairwise_phasing(input.seq = s,
                             input.twopt = tpt,
                             thresh.LOD.ph = 28,
                             max.conf.btnk.p1 = 20)
s.phased
err <- c(0.0, 0.01, 0.05)
L <- vector("list", length(err))
tm <- numeric(length(err))
for(i in 1:length(err))
  tm[i] <- system.time(L[[i]] <- hmm_map_reconstruction(input.seq = s.phased,
                                                        verbose = TRUE,
                                                        error = err[i],
                                                        tol = 10e-4))[3]
plot(tm~err, xlab = "error rate", ylab = "time(s)", col = 2, type = "b")
d <- lapply(L, function(x) sapply(x$phases, function(x) cumsum(c(0,imf_h(x$rf)))))
r <- c(0, max(sapply(d, max)))
plot(0, type = "n", ylim = r, xlim = c(0,n.mrk))
abline(h = map.length, lty = 2)
for(i in 1:length(err)){
  ll <- sapply(L[[i]]$phases, function(x) x$loglike)
  best <- which.max(ll)
  for(j in 1:ncol(d[[i]])){
    z1<-1; z2<-0.5
    if(j == best)
      z1 <- 3; z2 <- 1
      points(d[[i]][,j], type = "l", col = i+1, lwd = z1)
  }
}
legend("topleft", legend=err, col=1:length(err)+1, lty=1, cex=0.8)



i<-1
mrk.id <- rownames(s.phased$phases[[1]]$p1)
g <- s.phased$data$geno.dose[mrk.id, ]
id <- which(s.phased$data$ploidy.p2 == s.phased$data$dosage.p2[mrk.id])
g[id, ] <- g[id, ] - s.phased$data$ploidy.p2/2

bla<-mappoly2:::visit_states_biallelic_single(s.phased$phases[[i]]$p1, g, 0.01)


###### TEST #####
#i<-which.max(sapply(L[[1]]$phases, function(x) x$loglike))
G <- s.phased$data$geno.dose[mrk.id, ]
z1 <- mappoly2:::est_hmm_map_biallelic_single(s.phased$phases[[i]]$p1,
                                              G,
                                              rf = rep(0.01, (nrow(G)-1)),
                                              err = 0.05,
                                              verbose = T,
                                              detailed_verbose = F,
                                              tol = 10e-4,
                                              ret_H0 = F)
imf_h(z1[[2]])

#######
Z <- vector("list", length(s.phased$phases))
for(i in 1:length(s.phased$phases)){
  mrk.id <- rownames(s.phased$phases[[i]]$p1)
  G <- s.phased$data$geno.dose[mrk.id, ]
  pedigree <- matrix(rep(c(1,
                           2,
                           s.phased$data$ploidy.p1,
                           s.phased$data$ploidy.p2, 1),
                         s.phased$data$n.ind),
                     nrow = s.phased$data$n.ind,
                     byrow = TRUE)
  PH = list(s.phased$phases[[i]]$p1,
            s.phased$phases[[i]]$p2)
  Z[[i]] <- mappoly2:::est_hmm_map_biallelic2(PH, G, pedigree,
                                              rf = rep(0.01, (nrow(G)-1)),
                                              err = 0.0,
                                              verbose = T,
                                              detailed_verbose = F,
                                              tol = 10e-4,
                                              ret_H0 = F)
}

(best <- which.max(sapply(Z, function(x) x[[1]])))
imf_h(Z[[2]][[2]])
plot(c(0,cumsum(imf_h(Z[[best]][[2]]))))










### LOG algorithm ####
system.time(z2 <- mappoly2:::est_hmm_map_biallelic2(PH,
                                                    G, pedigree,
                                                    rf = rep(0.01, (nrow(G)-1)),
                                                    err = 0.05,
                                                    verbose = T,
                                                    detailed_verbose = F,
                                                    tol = 10e-4,
                                                    ret_H0 = F))
points(cumsum(imf_h(z1[[2]])), col = 4)
points(cumsum(imf_h(z2[[2]])), col = 2, pch = 20)

points(cumsum(imf_h(z1[[2]])), col = 3, pch = 4)
points(cumsum(imf_h(z2[[2]])), col = 5, pch = 5)


z1[[1]]
z2[[1]]




