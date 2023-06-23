rm(list = ls())
require(mappoly2)
source("~/repos/official_repos/misc/simulation.R")
ploidy.p1 = 2
ploidy.p2 = 4
n.mrk <- 100
ph<-test_simulate(ploidy.p1 = ploidy.p1,
                  ploidy.p2 = ploidy.p2,
                  fpath = "~/repos/official_repos/misc/fake_triploid.csv",
                  n.mrk = n.mrk,
                  n.ind = 200,
                  map.length =100,
                  miss.perc = 0,
                  n.chrom = 1,
                  random = FALSE,
                  seed = 279865)
dat <- read_geno_csv(file.in = "~/repos/official_repos/misc/fake_triploid.csv",
                     ploidy.p1 = ploidy.p1, ploidy.p2 = ploidy.p2)
s <- make_sequence(dat, "all", info.parent = "both")
tpt <- est_pairwise_rf(s)
s.phased <- pairwise_phasing(input.seq = s,
                             input.twopt = tpt,
                             thresh.LOD.ph = 15,
                             max.conf.btnk.p1 = 10)
s.phased
s.map <- hmm_map_reconstruction(input.seq = s.phased,
                                verbose = TRUE,
                                error = 0.01,
                                tol = 10e-4)
(ll <- sapply(s.map$phases, function(x) x$loglike))
which.max(ll)
x <- ll - max(ll)
names(x) <- seq_along(x)
x <- sort(x, decreasing = T)
plot(x, type = "n")
text(x, labels = names(x))
best <- which.max(ll)
d <- sapply(s.map$phases, function(x) cumsum(imf_h(x$rf)))
plot(d[,1], ylim = range(d), type = "l")
for(i in 2:ncol(d)){
  z<-1
  if(i == best)
    z <- 2
  points(d[,i], type = "l", col = z)
}
points(d[,best], type = "l", col = 3, lwd = 3)




##########
# (ll <- sapply(s.map$phases, function(x) x$loglike))

# for(i in 1:3){
#   w <- mappoly2:::vs_biallelic_error(PH = list(s.phased$phases[[i]]$p1,
#                                                   s.phased$phases[[i]]$p2),
#                                         G = g,
#                                         pedigree = pedigree,
#                                         err = 0.05)
#   print(w)
# }
#
#
# mrk.id <- rownames(s.phased$phases[[1]]$p1)
# g <- s.phased$data$geno.dose[mrk.id, ]
# pedigree <- matrix(rep(c(1,
#                          2,
#                          s.phased$data$ploidy.p1,
#                          s.phased$data$ploidy.p2, 1),
#                        s.phased$data$n.ind),
#                    nrow = s.phased$data$n.ind,
#                    byrow = TRUE)
# i<-1
# w <- mappoly2:::vs_biallelic_error(PH = list(s.phased$phases[[i]]$p1,
#                                              s.phased$phases[[i]]$p2),
#                                    G = g,
#                                    pedigree = pedigree,
#                                    err = 0.05)
#
#
#
# #l1 <- MAPs[[best]]$loglike
# l2 <- MAPs[[best]]$loglike
# p1.est <- mappoly::ph_matrix_to_list(PH[[best]]$p1)
# p1.sim <- mappoly::ph_matrix_to_list(ph[[1]]$p1[mrks,])
# mappoly::compare_haplotypes(dat$ploidy.p1, p1.est, p1.sim)
# p2.est <- mappoly::ph_matrix_to_list(PH[[best]]$p2)
# p2.sim <- mappoly::ph_matrix_to_list(ph[[1]]$p2[mrks,])
# mappoly::compare_haplotypes(dat$ploidy.p2, p2.est, p2.sim)
#
