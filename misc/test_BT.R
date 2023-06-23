rm(list = ls())
require(mappoly2)
require(crayon)
dat <- read_geno_csv(file.in = "~/repos/official_repos/misc/BT.csv",
                     ploidy.p1 = 6, ploidy.p2 = 6, name.p1 = "Beauregard",
                     name.p2 = "Tanzania")
dat <- filter_missing(dat, type = "marker", filter.thres = 0.2, inter = F)
dat <- filter_missing(dat, type = "individual", filter.thres = 0.2, inter = F)
dat <- filter_individuals(dat)
s.all <- filter_segregation(dat)
plot(s.all)
s9 <- make_sequence(dat, "ch09", info.parent = "p1")
s9
tpt9 <- est_pairwise_rf(s9, ncpus = 8)
m <- rf_list_to_matrix(tpt9,
                       shared.alleles = TRUE,
                       thresh.LOD.ph = 3)
g9 <- get_genome_order(s9)
s9 <- make_sequence(g9)
plot(m, ord = s9, fact = 5)
mrk.id <- s9$mrk.names
dose.p1 <- s9$data$dosage.p1[mrk.id]
ploidy.p1 <- s9$data$ploidy.p1
S.p1 <- m$Sh.p1[mrk.id, mrk.id]
Ph.p1 <- mappoly2:::twopt_phasing_cpp(mrk.id,ploidy.p1,dose.p1,S.p1, max_conf_number = 1)
dose.p2 <- s9$data$dosage.p2[mrk.id]
ploidy.p2 <- s9$data$ploidy.p2
S.p2 <- m$Sh.p2[mrk.id, mrk.id]
Ph.p2 <- mappoly2:::twopt_phasing_cpp(mrk.id,ploidy.p2,dose.p2,S.p2, max_conf_number = 2)
mrks <- intersect(Ph.p1$marker_names, Ph.p2$marker_names)
for(i in 1:length(Ph.p1$phase_configs))
  rownames(Ph.p1$phase_configs[[i]]) <- Ph.p1$marker_names
for(i in 1:length(Ph.p2$phase_configs))
  rownames(Ph.p2$phase_configs[[i]]) <- Ph.p2$marker_names
s9 <- make_sequence(dat, "ch09", info.parent = "p1")
mrks <- mrks[1:40]
s9 <- make_sequence(s9, mrks)
cte <- 1
PH <- MAPs <- vector("list", length(Ph.p1$phase_configs) * length(Ph.p2$phase_configs))
length(PH)
for(i in 1:length(Ph.p1$phase_configs)){
  for(j in 1:length(Ph.p2$phase_configs)){
    PH[[cte]] <- list(p1 = Ph.p1$phase_configs[[i]][mrks,],
                      p2 = Ph.p2$phase_configs[[j]][mrks,])
    MAPs[[cte]] <- hmm_map_reconstruction(s9,
                                          PH[[cte]],
                                          verbose = TRUE,
                                          tol = 10e-4)
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




