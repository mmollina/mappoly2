###### OLD ####
rm(list = ls())
require(mappoly2)
source("~/repos/official_repos/misc/simulation.R")
ploidy.p1 = 4
ploidy.p2 = 4
n.mrk <- 500
ph<-test_simulate(ploidy.p1 = ploidy.p1,
                  ploidy.p2 = ploidy.p2,
                  fpath = "~/repos/official_repos/misc/fake_triploid.csv",
                  n.mrk = n.mrk,
                  n.ind = 200,
                  map.length =100,
                  miss.perc = 0,
                  n.chrom = 1,
                  random = FALSE,
                  seed = 2986876)
dat <- read_geno_csv(file.in = "~/repos/official_repos/misc/fake_triploid.csv",
                     ploidy.p1 = ploidy.p1, ploidy.p2 = ploidy.p2)
s <- make_sequence(dat, "all", info.parent = "both")
tpt <- est_pairwise_rf(s, ncpus = 7)
s.phased <- pairwise_phasing(input.seq = s,
                             input.twopt = tpt,
                             thresh.LOD.ph = 5,
                             max.conf.btnk.p1 = 10)
s.phased

#### Working with error = 0.0 ########
s.map <- hmm_map_reconstruction(input.seq = s.phased,
                                verbose = TRUE,
                                error = 0.00,
                                tol = 10e-4)
best <- which.max(sapply(s.map$phases, function(x) x$loglike))
plot(cumsum(imf_h(s.map$phases[[best]]$rf)))

s.map <- hmm_map_reconstruction(input.seq = s.phased,
                                verbose = TRUE,
                                error = 0.01,
                                tol = 10e-4)
best <- which.max(sapply(s.map$phases, function(x) x$loglike))
points(cumsum(imf_h(s.map$phases[[best]]$rf)), col = 7, pch = 6)



