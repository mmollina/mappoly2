###### OLD ####
rm(list = ls())
require(mappoly2)
source("~/repos/official_repos/misc/simulation.R")
ploidy.p1 = 2
ploidy.p2 = 4
n.mrk <- 300
ph<-test_simulate(ploidy.p1 = ploidy.p1,
                  ploidy.p2 = ploidy.p2,
                  fpath = "~/repos/official_repos/misc/fake_triploid.csv",
                  n.mrk = n.mrk,
                  n.ind = 200,
                  map.length =100,
                  miss.perc = 10,
                  n.chrom = 1,
                  random = FALSE,
                  seed = 2986876)
####Read ####
dat <- read_geno_csv(file.in = "~/repos/official_repos/misc/fake_triploid.csv",
                     ploidy.p1 = ploidy.p1, ploidy.p2 = ploidy.p2)
dat
#### Select parent ####
s <- make_sequence(dat, "all", info.parent = "both")
#s <- make_sequence(dat, "all", info.parent = "p1")
#s <- make_sequence(dat, "all", info.parent = "p2")
#### Two-points ####
tpt <- est_pairwise_rf(s, ncpus = 7)

#### Phasing ####
s <- pairwise_phasing(input.seq = s,
                      input.twopt = tpt,
                      thresh.LOD.ph = 15,
                      max.conf.btnk.p1 = 10)
s
#### Mapping ####
s1 <- mapping(input.seq = s,
              verbose = TRUE,
              error = 0.00,# error = 0.0
              tol = 10e-4)
s1
s1 <- calc_haplotypes(s1)
l1 <- phase_remaining(s1,
                      tpt,
                      thresh.LOD.ph = 15,
                      max.phases = 10)
thresh.LOD.ph = 10
L<-l1[sapply(l1, function(x) length(x$loglike)==1 || x$loglike[2] > thresh.LOD.ph)]
names(L)











