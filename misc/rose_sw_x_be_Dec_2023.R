rm(list = ls())
require(mappoly2)
setwd("~/repos/official_repos/mappoly2/")
####Reading Data####
dat <- read_geno_csv(file.in = "misc/SWxBE.csv",
                     ploidy.p1 = 4,
                     ploidy.p2 = 4)
dat <- filter_data(dat, mrk.thresh = 0.2, ind.thresh = 0.19)
plot(dat, type = "density")
dat <- filter_individuals(dat)
#$209 seconds
system.time(dat <- pairwise_rf(dat,
                               mrk.scope = "per.chrom",
                               ncpus = 8))
plot(dat)
dat <- rf_filter(dat, probs = c(0.025, 0.975))
dat
g <- group(dat, expected.groups = 15, comp.mat = TRUE, inter = FALSE)
g
s <- make_sequence(g)
s <- order_sequence(s, type = "mds")
mappoly2:::plot_rf_matrix(s, type = "mds", fact = 5)
s <- rf_filter(s, type = "mds", probs = c(0.05, 1), diag.markers = 50)
mappoly2:::plot_rf_matrix(s, type = "mds", fact = 5)
##### Just parent 1 ####
s <- pairwise_phasing(s, type = "mds", parent = "p1", verbose = T)
s
s1 <- mapping(s1, parent = "p1", error = 0.0, tol = 10e-3, ncpus = 8)
map_summary(s1)
plot_map_list(s1)
