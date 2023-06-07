###### OLD ####
rm(list = ls())
require(mappoly2)
source("~/repos/official_repos/misc/simulation.R")
ploidy.p1 = 4
ploidy.p2 = 4
n.mrk <- 200
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
s <- make_sequence(dat, "all", info.parent = "p1")
tpt <- est_pairwise_rf(s)
s.phased <- pairwise_phasing(input.seq = s,
                             input.twopt = tpt,
                             thresh.LOD.ph = 5,
                             max.conf.btnk.p1 = 10)
s.phased
#### error = 0.0 ########
s.map1 <- hmm_map_reconstruction(input.seq = s.phased,
                                 verbose = TRUE,
                                 error = 0.00,
                                 tol = 10e-4)
best <- which.max(sapply(s.map1$phases, function(x) x$loglike))
plot(cumsum(imf_h(s.map1$phases[[best]]$rf)))
#### error = 0.05 ########
s.map2 <- hmm_map_reconstruction(input.seq = s.phased,
                                 verbose = TRUE,
                                 error = 0.05,
                                 tol = 10e-4)
best <- which.max(sapply(s.map2$phases, function(x) x$loglike))
points(cumsum(imf_h(s.map2$phases[[best]]$rf)), col = 7, pch = 6)
#########################################
mrk.id <- rownames(s.map1$phases[[best]]$p1)
g <- s.map1$data$geno.dose[mrk.id, ]
pedigree <- matrix(rep(c(1,
                         2,
                         s.map1$data$ploidy.p1,
                         s.map1$data$ploidy.p2, 1),
                       s.map1$data$n.ind),
                   nrow = s.map1$data$n.ind,
                   byrow = TRUE)
PH = list(s.map1$phases[[best]]$p1,
          s.map1$phases[[best]]$p2)
w1 <- mappoly2:::calc_genoprob_biallelic(PH = PH,
                                         G = g,
                                         pedigree = pedigree,
                                         rf = s.map1$phases[[best]]$rf,
                                         err = 0.0)
w2 <- mappoly2:::calc_genoprob_biallelic(PH = PH,
                                         G = g,
                                         pedigree = pedigree,
                                         rf = s.map2$phases[[best]]$rf,
                                         err = 0.05)
i<-4
layout(matrix(1:2,2,1))
image(w1[[i]])
image(w2[[i]])
#############################################
mrk.id <- rownames(s.map1$phases[[best]]$p1)
g <- s.map1$data$geno.dose[mrk.id, ]
pedigree <- matrix(rep(c(1,
                         2,
                         s.map1$data$ploidy.p1,
                         s.map1$data$ploidy.p2, 1),
                       s.map1$data$n.ind),
                   nrow = s.map1$data$n.ind,
                   byrow = TRUE)
PH = s.map1$phases[[best]]$p1
id <- which(s.map1$data$ploidy.p2 == s.map1$data$dosage.p2[mrk.id])
g[id, ] <- g[id, ] - s.map1$data$ploidy.p2/2
w1 <- mappoly2:::calc_genoprob_biallelic_single(PH = PH,
                                                G = g,
                                                rf = s.map1$phases[[best]]$rf,
                                                err = 0.0)
w2 <- mappoly2:::calc_genoprob_biallelic_single(PH = PH,
                                                G = g,
                                                rf = s.map2$phases[[best]]$rf,
                                                err = 0.05)
i<-4
layout(matrix(1:2,2,1))
image(w1[[i]])
image(w2[[i]])

