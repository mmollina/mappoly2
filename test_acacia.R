rm(list = ls())
require(mappoly2)
P1 <- "AC46"
P2 <- "AC78"
x <- read_geno_csv(file.in = paste0("~/repos/collaborations/acacia/full_sib_data/mappoly2_", P1, "-", P2,".csv"),
                     ploidy.p1 = 2,
                     ploidy.p2 = 2,
                     name.p1 = P1,
                     name.p2 = P2)
x <- filter_data(x, mrk.thresh = .08, ind.thresh = .7)
x <- set_initial_sequence(x, "ch1")
system.time(x <- pairwise_rf(x, ncpus = 7))
x <- init_rf_filter(x, 5, 5, .1, c(0.01, 0.99))
plot(x, what = "pairwise", fact = 1, ord = x$initial.screened.rf$mrk.names)
x <- set_working_sequence(x, ch = list(1))
x <- genome_order_per_group(x, gr = 1)
plot(x, what = "pairwise", ord = rownames(x$working.sequences[[i]]$order$genome$info))
x <- mds_per_group(x, gr = 1)
plot(x$working.sequences[[i]]$order$mds$info$locimap$confplotno)
plot(x, what = "pairwise", ord = x$working.sequences[[i]]$order$mds$info$locimap$locus)
x <- pairwise_phasing_per_group(x,
                                gr = 1,
                                type = "genome",
                                thresh.LOD.ph = 20,
                                thresh.LOD.rf = 20)
x <- mapping_per_group(x, type = "genome", gr = 1, error = 0.05)
plot_map(x, gr = 1, type = "genome")


