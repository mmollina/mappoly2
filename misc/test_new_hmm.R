###### OLD ####
rm(list = ls())
require(mappoly2)
source("~/repos/official_repos/misc/simulation.R")
ploidy.p1 = 4
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
                             thresh.LOD.ph = 5,
                             max.conf.btnk.p1 = 10)
s.phased
s.map <- hmm_map_reconstruction(input.seq = s.phased,
                                verbose = TRUE,
                                tol = 10e-4)
plot(cumsum(imf_h(s.map$phases[[1]]$rf)))

###### TEST #####
mrk.id <- rownames(s.phased$phases[[1]]$p1)
G <- s.phased$data$geno.dose[mrk.id, ]
pedigree <- matrix(rep(c(1,
                         2,
                         s.phased$data$ploidy.p1,
                         s.phased$data$ploidy.p2, 1),
                       s.phased$data$n.ind),
                   nrow = s.phased$data$n.ind,
                   byrow = TRUE)
i<-1
PH = list(s.phased$phases[[i]]$p1,
          s.phased$phases[[i]]$p2)
w <- mappoly2:::vs_biallelic_error2(PH, G, pedigree, 0.01)
w$emit[[1]][[1]]
z <- mappoly2:::est_hmm_map_biallelic2(PH, G, pedigree,
                                       rf = rep(0.01, (nrow(G)-1)),
                                       verbose = F,
                                       detailed_verbose = F,
                                       tol = 10e-3,
                                       ret_H0 = F)
