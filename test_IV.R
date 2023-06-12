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
                      max.conf.btnk.p1 = 20)
s
#### Mapping ####
s1 <- mapping(input.seq = s,
              verbose = TRUE,
              error = 0.00,# error = 0.0
              tol = 10e-4)

plot(cumsum(imf_h(s1$phases[[1]]$rf)))

x<- NULL
for(i in seq(30,0,-5)){
  m <- rf_list_to_matrix(tpt,
                         thresh.LOD.ph = i,
                         shared.alleles = TRUE)
  mrk.pos <- rownames(s1$phases[[1]]$p2)
  mrk.id <- setdiff(s1$mrk.names, mrk.pos)
  InitPh <- s1$phases[[1]]$p2
  dose.vec <- s1$data$dosage.p2[mrk.id]
  S <- m$Sh.p2[mrk.id, mrk.pos]
  verbose <- TRUE
  L <- mappoly2:::phasing_one(mrk.id, dose.vec, S, InitPh, verbose)
  x <- cbind(x, sort(sapply(L, nrow)))
}
image(x)
