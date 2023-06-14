###### OLD ####
rm(list = ls())
require(mappoly2)
source("~/repos/official_repos/misc/simulation.R")
ploidy.p1 = 4
ploidy.p2 = 4
n.mrk <- 400
ph<-test_simulate(ploidy.p1 = ploidy.p1,
                  ploidy.p2 = ploidy.p2,
                  fpath = "~/repos/official_repos/misc/fake_triploid.csv",
                  n.mrk = n.mrk,
                  n.ind = 300,
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
tpt <- est_pairwise_rf(s)

#### Phasing ####
s <- pairwise_phasing(input.seq = s,
                      input.twopt = tpt,
                      thresh.LOD.ph = 5,
                      max.conf.btnk.p1 = 20)
s
#### Mapping ####
s1 <- mapping(input.seq = s,
              verbose = TRUE,
              error = 0.00,# error = 0.0
              tol = 10e-4)
#s2 <- mapping(input.seq = s,
#              verbose = TRUE,
#              error = 0.05,# error = 0.05
#              tol = 10e-4)

plot(cumsum(imf_h(s1$phases[[1]]$rf)))
#points(cumsum(imf_h(s2$phases[[1]]$rf)), col = 7, pch = 6)

#### Haplotypes ####
s1 <- calc_haplotypes(s1)
#s2 <- calc_haplotypes(s2)

image((as.matrix(s1$phases[[1]]$haploprob[1:500, -c(1:2)])), main = s1$phases[[1]]$error)
#image((as.matrix(s2$phases[[1]]$haploprob[1:16, -c(1:2)])), main = s2$phases[[1]]$error)

#x<-s1$phases[[1]]$haploprob[1:6, 6]
#x<-c(0.2,0.8,0,1,1,0)
#mappoly2:::homologprob_to_hmmstates(h = x,
#                                    ploidy1 = 2, ploidy2 = 4)
#s1$phases[[1]]$haploprob[1:12,1:6]

z<-mappoly2:::homologprob_to_hmmstates(H = s1$phases[[1]]$haploprob, rep(6,300), rep(6,300), k = 4)








