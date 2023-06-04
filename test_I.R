###### OLD ####
rm(list = ls())
require(mappoly2)
source("~/repos/official_repos/misc/simulation.R")
ploidy.p1 = 4
ploidy.p2 = 4
n.mrk <- 10
ph<-test_simulate(ploidy.p1 = ploidy.p1,
                  ploidy.p2 = ploidy.p2,
                  fpath = "~/repos/official_repos/misc/fake_triploid.csv",
                  n.mrk = n.mrk,
                  n.ind = 200,
                  map.length =10,
                  miss.perc = 0,
                  n.chrom = 1,
                  random = FALSE,
                  seed = 2986)
dat <- read_geno_csv(file.in = "~/repos/official_repos/misc/fake_triploid.csv",
                     ploidy.p1 = ploidy.p1, ploidy.p2 = ploidy.p2)
s <- make_sequence(dat, "all", info.parent = "both")
tpt <- est_pairwise_rf(s)
s.phased <- pairwise_phasing(input.seq = s,
                             input.twopt = tpt,
                             thresh.LOD.ph = 3,
                             max.conf.btnk.p1 = 2)
s.phased
s.map <- hmm_map_reconstruction(input.seq = s.phased,
                                verbose = TRUE,
                                tol = 10e-4)
plot(cumsum(imf_h(s.map$phases[[1]]$rf)))

###### TEST #####
i<-1
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
w <- mappoly2:::vs_biallelic_error(PH, G, pedigree, 0.0)
w$emit[[1]][[1]]
w$states[[1]][[1]]
system.time(z1 <- mappoly2:::est_hmm_map_biallelic(PH, G, pedigree,
                                       rf = rep(0.01, (nrow(G)-1)),
                                       verbose = T,
                                       detailed_verbose = F,
                                       tol = 10e-4,
                                       ret_H0 = F))
system.time(z2 <- mappoly2:::est_hmm_map_biallelic2(PH, G, pedigree,
                                                   rf = rep(0.01, (nrow(G)-1)),
                                                   verbose = T,
                                                   detailed_verbose = F,
                                                   tol = 10e-4,
                                                   ret_H0 = F))

round(imf_h(z1[[2]]), 3)
round(imf_h(z2[[2]]), 3)
z1[[1]]
z2[[1]]




