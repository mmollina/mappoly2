###### OLD ####
rm(list = ls())
require(mappoly2)
source("~/repos/official_repos/misc/simulation.R")
ploidy.p1 = 2
ploidy.p2 = 4
n.mrk <- 700
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
s.phased <- pairwise_phasing(input.seq = s,
                             input.twopt = tpt,
                             thresh.LOD.ph = 5,
                             max.conf.btnk.p1 = 10)
s.phased
#### Mapping ####
s.map1 <- mapping(input.seq = s.phased,
                                verbose = TRUE,
                                error = 0.00,# error = 0.0
                                tol = 10e-4)
s.map2 <- mapping(input.seq = s.phased,
                                  verbose = TRUE,
                                  error = 0.05,# error = 0.05
                                  tol = 10e-4)

plot(cumsum(imf_h(s.map1$phases[[1]]$rf)))
points(cumsum(imf_h(s.map2$phases[[1]]$rf)), col = 7, pch = 6)
#### Haplotypes ####
s.haplo1 <- calc_haplotypes(s.map1)
s.haplo2 <- calc_haplotypes(s.map2)

image((as.matrix(s.haplo1$phases[[1]]$haploprob[1:16, -c(1:2)])), main = s.haplo1$phases[[1]]$error)
image((as.matrix(s.haplo2$phases[[1]]$haploprob[1:16, -c(1:2)])), main = s.haplo2$phases[[1]]$error)


