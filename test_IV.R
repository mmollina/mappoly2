###### OLD ####
rm(list = ls())
require(mappoly2)
source("~/repos/official_repos/misc/simulation.R")
ploidy.p1 = 4
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
                  seed = 26876)
####Read ####
dat <- read_geno_csv(file.in = "~/repos/official_repos/misc/fake_triploid.csv",
                     ploidy.p1 = ploidy.p1, ploidy.p2 = ploidy.p2)
dat
#### Select parent ####
s <- make_sequence(dat, 1:300, info.parent = "both")
#s <- make_sequence(dat, "all", info.parent = "p1")
#s <- make_sequence(dat, "all", info.parent = "p2")
#### Two-points ####
tpt <- est_pairwise_rf(s, ncpus = 7)

#### Phasing ####
s <- pairwise_phasing(input.seq = s,
                      input.twopt = tpt,
                      thresh.LOD.ph = 5,
                      max.conf.btnk.p1 = 5)
s
#### Mapping ####
s1 <- mapping(input.seq = s,
              verbose = TRUE,
              error = 0.00,# error = 0.0
              tol = 10e-4)
plot_map(s1)
s1 <- calc_haplotypes(s1)
s1 <- augment_phased_map(input.seq = s1,
                         input.twopt = tpt,
                         thresh.LOD.ph = 15,
                         thresh.LOD.rf = 3,
                         thresh.rf = 0.5,
                         max.phases = 10,
                         thresh.LOD.ph.to.insert = 10,
                         tol = 10e-4,
                         verbose = TRUE)
s1 <- mapping(input.seq = s1,
              verbose = TRUE,
              error = 0.00,# error = 0.0
              tol = 10e-4)
s1 <- calc_haplotypes(s1)
plot_map(s1, left.lim = 10, right.lim = 40, mrk.names = T)








