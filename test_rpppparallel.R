rm(list = ls())
require(mappoly2)
setwd("/Users/mmollin/repos/official_repos/mappoly2")
source("misc/simulation.R")
ploidy.p1 = 4
ploidy.p2 = 4
n.mrk <- 100
ph<-test_simulate(ploidy.p1 = ploidy.p1,
                  ploidy.p2 = ploidy.p2,
                  fpath = "misc/fake_triploid.csv",
                  n.mrk = n.mrk,
                  n.ind = 200,
                  map.length = 100,
                  miss.perc = 0,
                  n.chrom = 3,
                  random = FALSE,
                  seed = 43598)
dat <- read_geno_csv(file.in = "misc/fake_triploid.csv",
                     ploidy.p1 = ploidy.p1,
                     ploidy.p2 = ploidy.p2,
                     name.p1 = "parent_1",
                     name.p2 = "parent_2")
s <- make_sequence(dat, "all")
s <- pairwise_rf(input.seq = s, ncpus = 2)
plot(s, "rf")
