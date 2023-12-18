rm(list = ls())
require(mappoly2)
setwd("~/repos/official_repos/mappoly2/")
####Reading Data####
dat <- read_geno_csv(file.in = "misc/I195xJ432.csv",
                     ploidy.p1 = 4,
                     ploidy.p2 = 4)
print(dat)
plot(dat)

#### Initial QA/QC ####
dat <- filter_data(dat, mrk.thresh = 0.2, ind.thresh = 0.1)
plot(dat)
dat
plot(dat, type = "density")
dat <- filter_individuals(dat, ind.to.remove = "F1.85.133")
#dat <- filter_individuals(dat, type = "PCA")
dat
#### Pairwise rf ####
## 53 seconds
system.time(dat <- pairwise_rf(dat, mrk.scope = "all", ncpus = 8))
dat
plot(dat)
#### RF-based filter ####
dat <- rf_filter(dat, probs = c(0.025, 0.975))
dat
#plot(dat)

#### Grouping ####
g <- group(x = dat, expected.groups = 8, comp.mat = TRUE, inter = FALSE)
g
#### Sequence ####
s <- make_sequence(g,
                   ch = list(1, 2, 3, 4, 5, 6, 7, 8),
                   lg = list(1, 3, 2, 4, 7, 6, 5, 8))
s
# g <- group(x = dat, expected.groups = 15, comp.mat = TRUE, inter = FALSE)
# g
# s <- make_sequence(g,
#                    ch = list(1, 2, 3, 4, 5, 6, 7, 8),
#                    lg = list(1,
#                              c(3,7),
#                              c(2,5),
#                              c(4,11,15),
#                              10,
#                              c(8,14),
#                              c(6,9),
#                              c(12,13)))
#s
#s <- make_sequence(g, lg = list(1,3,2,4,6,7,5,8))

#### Order ######
s <- order_sequence(s, type = "mds")
s
s <- order_sequence(s, type = "genome")
print(s, type = "genome")
#mappoly2:::plot_rf_matrix(s, type = "mds")
#mappoly2:::plot_rf_matrix(s, type = "genome")

#### RF-based filter per groups ####
s <- rf_filter(s, type = "mds", probs = c(0.05, 1))
s <- rf_filter(s, type = "genome", probs = c(0.05, 1), diag.markers = 50)
#mappoly2:::plot_rf_matrix(s, type = "mds")
#mappoly2:::plot_rf_matrix(s, type = "genome")

#### Pairwise Phasing ####
s <- pairwise_phasing(s, type = "mds", parent = "p1")
s
s <- pairwise_phasing(s, type = "mds", parent = "p2")
s
s <- pairwise_phasing(s, type = "mds", parent = "p1p2")
s
s <- pairwise_phasing(s, type = "genome", parent = "p1")
s <- pairwise_phasing(s, type = "genome", parent = "p2")
s <- pairwise_phasing(s, type = "genome", parent = "p1p2")
print(s, type = "genome")
#### Mapping MDS####
system.time(s <- mapping(s,
                         type = "mds",
                         parent = "p1",
                         tol = 10e-4,
                         error = 0.0,
                         ncpus = 1))
system.time(s <- mapping(s,
                         type = "mds",
                         parent = "p2",
                         tol = 10e-4,
                         error = 0.0,
                         ncpus = 1))
system.time(s <- mapping(s,
                         type = "mds",
                         parent = "p1p2",
                         tol = 10e-4,
                         error = 0.0,
                         ncpus = 8))
plot_map_list(s, parent = "p1")
plot_map_list(s, parent = "p2")
plot_map_list(s, parent = "p1p2")
plot_map(s, parent = "p1", lg = 1)
plot_map(s, parent = "p2", lg = 1)
plot_map(s, parent = "p1p2", lg = 1)
mappoly2:::map_summary(s, type = "mds", parent = "p1")
mappoly2:::map_summary(s, type = "mds", parent = "p2")
mappoly2:::map_summary(s, type = "mds", parent = "p1p2")
#### Mapping Genome####
system.time(s <- mapping(s,
                         type = "genome",
                         parent = "p1",
                         tol = 10e-4,
                         error = 0.0,
                         ncpus = 1))
system.time(s <- mapping(s,
                         type = "genome",
                         parent = "p2",
                         tol = 10e-4,
                         error = 0.0,
                         ncpus = 1))
plot_map_list(s, type = "genome",parent = "p1")
plot_map_list(s, type = "genome",parent = "p2")
plot_map(s, type = "genome", parent = "p1", lg = 6)
plot_map(s, type = "genome", parent = "p2", lg = 1)
mappoly2:::map_summary(s, type = "genome", parent = "p1")
print(s, type = "genome")
#### Haplotypes ####
s <- calc_haplotypes(s, type = "mds", ncpus = 1)
s <- calc_haplotypes(s, type = "genome", ncpus = 1, parent = "p1")
s <- calc_haplotypes(s, type = "genome", ncpus = 1, parent = "p2")
print(s, type = "genome")
#### Merging P1 and P2 (will overwrite p1p2 slot) ####
s <- merge_single_parent_maps(s, type = "genome", ncpus = 8)
print(s, "genome")
plot_map_list(s, type = "genome",parent = "p1p2", col = viridis::turbo(8, begin = .2, end = .8))
s <- mapping(s, type = "genome", error = 0.05, ncpus = 8)
print(s, "genome")
plot_map_list(s, type = "genome",parent = "p1p2", col = viridis::turbo(8, begin = .2, end = .8))
s <- calc_haplotypes(s, type = "genome", ncpus = 1, parent = "p1p2")
print(s, "genome")
map_summary(s, "genome")
plot_map(s, lg = 1, type = "genome")
#### Augment maps ####








plot_genome_vs_map(s)
plot_genome_vs_map(s, same.ch.lg = T, alpha = 1, size = 2)
plot_genome_vs_map(s, same.ch.lg = T, alpha = 1, size = 2, type = "genome")

x1<-s$maps$lg1$mds$phase[[1]]$haploprob
x1[1:10, 1:10]
image(t(as.matrix(x1[1:8,-c(1:3)])))

x2<-s$maps$lg1$genome$phase[[1]]$haploprob
round(x2[1:10, 1:10],2)
image(t(as.matrix(x2[1:8,-c(1:3)])))

plot_genome_vs_map(s)
s <- rev_map(s, 2:8)
plot_genome_vs_map(s)

##### Just parent 1 ####
s1 <- pairwise_phasing(s, type = "mds", parent = "p1", verbose = F)
s1
s1 <- mapping(s1, parent = "p1", error = 0.05)
map_summary(s1)
plot_map_list(s1)


plot_map(s, lg = 1, type = "mds", mrk.names = TRUE)
plot_map(s, lg = 2, type = "mds", mrk.names = TRUE)
plot_map(s, lg = 3, type = "mds", mrk.names = TRUE)
plot_map(s, lg = 4, type = "mds", mrk.names = TRUE)
plot_map(s, lg = 5, type = "mds", mrk.names = TRUE)
plot_map(s, lg = 6, type = "mds", mrk.names = TRUE)
plot_map(s, lg = 7, type = "mds", mrk.names = TRUE)
plot_map(s, lg = 8, type = "mds", mrk.names = TRUE)


s.ch1.all <- mapping(input.seq = s.ch1.all,
                     verbose = TRUE,
                     error = 0.00,
                     tol = 10e-4)
plot_map(s.ch1.all, mrk.names = TRUE)
s.ch1.all <- calc_haplotypes(s.ch1.all)
s.ch1.all <- augment_phased_map(input.seq = s.ch1.all,
                                input.twopt = tpt)
s.ch1.all <- mapping(input.seq = s.ch1.all,
                     verbose = TRUE,
                     error = 0.00,
                     tol = 10e-4)
plot_map(s.ch1.all)
s.ch1.all <- calc_haplotypes(s.ch1.all)



### Begin
mrk <- "Ch_1_M_1"
drop.seq <- drop_marker(input.seq = s.ch1.all, mrk, reestimate.map = T)
plot_map(drop.seq, 0, 10, mrk.names = T, xlim = c(-2, 11), cex = 1.5)
drop.seq<-calc_haplotypes(drop.seq)
drop.seq <- add_marker(drop.seq, mrk, thresh.rf.to.add = .1)
plot_map(s.ch1.all, 0, 10, mrk.names = T, xlim = c(-2, 11), cex = 1.5)
x11()
plot_map(drop.seq, 0, 10, mrk.names = T, xlim = c(-2, 11), cex = 1.5)
drop.seq <- mapping(drop.seq, verbose = T)
plot_map(drop.seq, 0, 10, mrk.names = T, xlim = c(-2, 11), cex = 1.5)

### Middle
mrk <- "Ch_1_M_74"
drop.seq <- drop_marker(input.seq = s.ch1.all, mrk, reestimate.map = T)
plot_map(drop.seq, 60, 80, mrk.names = T, xlim = c(58, 81), cex = 1.5)
drop.seq<-calc_haplotypes(drop.seq)
drop.seq <- add_marker(drop.seq, mrk)
plot_map(drop.seq, 60, 80, mrk.names = T, xlim = c(58, 81), cex = 1.5)
drop.seq <- mapping(drop.seq, verbose = T)
plot_map(drop.seq, 60, 80, mrk.names = T, xlim = c(58, 81), cex = 1.5)


#### End
plot_map(s.ch1.all, 90, 100, mrk.names = T, xlim = c(89, 101), cex = 1.5)
mrk <- "Ch_1_M_100"
drop.seq <- drop_marker(input.seq = s.ch1.all, mrk, reestimate.map = T)
plot_map(drop.seq, 90, 100, mrk.names = T, xlim = c(89, 101), cex = 1.5)
drop.seq<-calc_haplotypes(drop.seq)
drop.seq <- add_marker(drop.seq, mrk, thresh.rf.to.add = 0.1, tol = 10e-6)
plot_map(drop.seq, 90, 100, mrk.names = T, xlim = c(89, 103), cex = 1.5)
drop.seq <- mapping(drop.seq, verbose = T)
plot_map(drop.seq, 90, 100, mrk.names = T, xlim = c(89, 101), cex = 1.5)























