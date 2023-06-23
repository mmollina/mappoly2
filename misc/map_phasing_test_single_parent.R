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
                  n.ind = 1000,
                  map.length = 100,
                  miss.perc = 0,
                  n.chrom = 1,
                  random = FALSE,
                  seed = 274985)
dat <- read_geno_csv(file.in = "~/repos/official_repos/misc/fake_triploid.csv",
                     ploidy.p1 = ploidy.p1, ploidy.p2 = ploidy.p2)
s1 <- make_sequence(dat, "all",  info.parent = "p1")
mrk.id <- s1$mrk.names
tpt <- est_pairwise_rf(s1)
m <- rf_list_to_matrix(tpt, shared.alleles = T, thresh.LOD.ph = 1)
image(m$Sh.p1)
dose.p1 <- s1$data$dosage.p1[mrk.id]
ploidy.p1 <- s1$data$ploidy.p1
S.p1 <- m$Sh.p1[mrk.id, mrk.id]
Ph.p1 <- mappoly2:::twopt_phasing_cpp(mrk.id,ploidy.p1,dose.p1,S.p1, max_conf_number = 10)
dose.p2 <- s1$data$dosage.p2[mrk.id]
ploidy.p2 <- s1$data$ploidy.p2
S.p2 <- m$Sh.p2[mrk.id, mrk.id]
Ph.p2 <- mappoly2:::twopt_phasing_cpp(mrk.id,ploidy.p2,dose.p2,S.p2, max_conf_number = 4)
mrks <- intersect(Ph.p1$marker_names, Ph.p2$marker_names)
for(i in 1:length(Ph.p1$phase_configs))
  rownames(Ph.p1$phase_configs[[i]]) <- Ph.p1$marker_names
for(i in 1:length(Ph.p2$phase_configs))
  rownames(Ph.p2$phase_configs[[i]]) <- Ph.p2$marker_names
s1 <- make_sequence(s1, mrks)
cte <- 1
PH <- MAPs <- vector("list", length(Ph.p1$phase_configs) * length(Ph.p2$phase_configs))
length(PH)
for(i in 1:length(Ph.p1$phase_configs)){
  for(j in 1:length(Ph.p2$phase_configs)){
    PH[[cte]] <- list(p1 = Ph.p1$phase_configs[[i]][mrks,],
                      p2 = Ph.p2$phase_configs[[j]][mrks,])
    MAPs[[cte]] <- hmm_map_reconstruction(s1,PH[[cte]],verbose = TRUE,tol = 10e-4)
    cte <- cte + 1
  }
}
best <- which.max(sapply(MAPs, function(x) x$loglike))
plot(cumsum(imf_h(MAPs[[best]]$rec.frac)))
p1.est <- mappoly::ph_matrix_to_list(PH[[best]]$p1)
p1.sim <- mappoly::ph_matrix_to_list(ph[[1]]$p1[mrks,])
mappoly::compare_haplotypes(dat$ploidy.p1, p1.est, p1.sim)
p2.est <- mappoly::ph_matrix_to_list(PH[[best]]$p2)
p2.sim <- mappoly::ph_matrix_to_list(ph[[1]]$p2[mrks,])
mappoly::compare_haplotypes(dat$ploidy.p2, p2.est, p2.sim)

