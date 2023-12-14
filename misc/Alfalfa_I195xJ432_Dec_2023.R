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
plot(dat)

#### Grouping ####
g <- group(dat, expected.groups = 8, comp.mat = TRUE, inter = FALSE)
g
#### Sequence ####
s <- make_sequence(g,
                   ch = list(1, 3, 2, 4, 7, 6, 5, 8),
                   lg = list(1, 2, 3, 4, 5, 6, 7, 8),
                   seq.order = c(1, 3, 2, 4, 7, 6, 5, 8))
#s
#s <- make_sequence(g, lg = list(1,3,2,4,6,7,5,8))

#### Order ######
s <- order_sequence(s, type = "mds")
s
s <- order_sequence(s, type = "genome")
s
mappoly2:::plot_rf_matrix(s, type = "mds")
mappoly2:::plot_rf_matrix(s, type = "genome")

#### RF-based filter per groups ####
s <- rf_filter(s, type = "mds", probs = c(0.05, 1))
s <- rf_filter(s, type = "genome", probs = c(0.05, 1), diag.markers = 50)
mappoly2:::plot_rf_matrix(s, type = "mds")
mappoly2:::plot_rf_matrix(s, type = "genome")

#### Pairwise Phasing ####
s <- pairwise_phasing(s, type = "mds")
s
s <- pairwise_phasing(s, type = "genome")
s
#### Mapping ####
system.time(s <- mapping(s,
              type = "mds",
              tol = 10e-3,
              error = 0.0,
              ncpus = 8))
system.time(s <- mapping(s,
             type = "mds",
             tol = 10e-3,
             error = 0.05,
             ncpus = 8))
s <- mapping(s,
             type = "genome",
             tol = 10e-3,
             error = 0.0,
             ncpus = 8)
s <- mapping(s,
             type = "genome",
             tol = 10e-3,
             error = 0.05,
             ncpus = 8)
mappoly2:::map_summary(s, type = "mds")
mappoly2:::map_summary(s, type = "genome")

mappoly2:::plot_map_list(s,horiz = FALSE, col = mappoly:::mp_pallet3(8))
mappoly2:::plot_map_list(s,horiz = TRUE, col = mappoly:::mp_pallet3(8))
mappoly2:::plot_map_list(s, type = "genome", horiz = FALSE, col = mappoly:::mp_pallet3(8))
mappoly2:::plot_map_list(s, type = "genome", horiz = TRUE, col = mappoly:::mp_pallet3(8))


plot_map(s, lg = 1, type = "mds", mrk.names = TRUE)
plot_map(s, lg = 2, type = "mds", mrk.names = TRUE)
plot_map(s, lg = 3, type = "mds", mrk.names = TRUE)
plot_map(s, lg = 4, type = "mds", mrk.names = TRUE)
plot_map(s, lg = 5, type = "mds", mrk.names = TRUE)
plot_map(s, lg = 6, type = "mds", mrk.names = TRUE)
plot_map(s, lg = 7, type = "mds", mrk.names = TRUE)
plot_map(s, lg = 8, type = "mds", mrk.names = TRUE)

plot_map(s, lg = 1, type = "genome", mrk.names = TRUE)
plot_map(s, lg = 2, type = "genome", mrk.names = TRUE)
plot_map(s, lg = 3, type = "genome", mrk.names = TRUE)
plot_map(s, lg = 4, type = "genome", mrk.names = TRUE)
plot_map(s, lg = 5, type = "genome", mrk.names = TRUE)
plot_map(s, lg = 6, type = "genome", mrk.names = TRUE)
plot_map(s, lg = 7, type = "genome", mrk.names = TRUE)
plot_map(s, lg = 8, type = "genome", mrk.names = TRUE)

#### Haplotypes ####
s <- calc_haplotypes(s, type = "mds", ncpus = 1)
s <- calc_haplotypes(s, type = "mds", ncpus = 1)
s <- calc_haplotypes(s, type = "genome", ncpus = 1)

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


#save(s, file = "alfalfa_map.rda")



#### Select parent 2####
s.ch1.p2 <- make_sequence(s.ch1.all, info.parent = "p2")
s.ch1.p2
# Phasing #
s.ch1.p2 <- pairwise_phasing(input.seq = s.ch1.p2,
                             input.twopt = tpt,
                             thresh.LOD.ph = 3,
                             max.conf.btnk.p1 = 10)
# Mapping #
s.ch1.p2 <- mapping(input.seq = s.ch1.p2,
                    verbose = TRUE,
                    error = 0.00,# error = 0.0
                    tol = 10e-4)
plot_map(s.ch1.p2, mrk.names = TRUE)
#### Merging maps from p1 and p2 ####
s.ch1.all <- merge_single_parent_maps(input.seq.all = s.ch1.all,
                                      input.seq.p1 = s.ch1.p1,
                                      input.seq.p2 = s.ch1.p2,
                                      input.twopt = tpt.ch1)
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
input.mat <- rf_list_to_matrix(tpt.ch1,
                               thresh.LOD.ph = 5,
                               thresh.LOD.rf = 5,
                               thresh.rf = .5,
                               shared.alleles = TRUE)
plot_map(s.ch1.all, mrk.names = T, xlim = c(-10,105))



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



























x <- set_initial_sequence(x, "all")
plot(x, what = "initiated")
x <- set_initial_sequence(x, c(1:100))
plot(x, what = "initiated")
x <- set_initial_sequence(x, sample(x$data$mrk.names, 100))
plot(x, what = "initiated")
x <- set_initial_sequence(x, c("ch3"))
plot(x, what = "initiated")

#### Compute RF ####
x <- set_initial_sequence(x, "all")
n.cpus <- parallel::detectCores() - 1
x <- pairwise_rf(x, ncpus = n.cpus)



plot(x, what = "pairwise", fact = 5)
plot(x, what = "pairwise", type = "lod", fact = 5)
#### Initial filter based in RF ####
x <- init_rf_filter(x, 5, 5, .1, c(0.01, 0.99))
plot(x, what = "pairwise")

#### Grouping ####
x <- group(x, expected.groups = 13)
plot(x, what = "group")
print(x)
#### Set Working Sequences ####
x <- set_working_sequence(x,
                          lg = list(c(1,4),
                                    c(2,3),
                                    c(5,6)),
                          ch = list(c(1,3),
                                    2,
                                    c(4,5)))

x <- set_working_sequence(x, ch = list(1,2,3,4,5))
#### Working sequence filter based in RF ####
x <- rf_filter_per_group(x, gr = 1, 5, 5, .1, c(0.01, 0.98))
x <- rf_filter_per_group(x, gr = 2, 5, 5, .1, c(0.01, 0.98))
x <- rf_filter_per_group(x, gr = 3, 5, 5, .1, c(0.01, 0.98))
x <- rf_filter_per_group(x, gr = 4, 5, 5, .1, c(0.04, 0.98))
x <- rf_filter_per_group(x, gr = 5, 5, 5, .1, c(0.01, 0.98))
#### Ordering ####
x <- mds_per_group(x, gr = 1)
x <- mds_per_group(x, gr = 2)
x <- mds_per_group(x, gr = 3)
x <- mds_per_group(x, gr = 4)
x <- mds_per_group(x, gr = 5)
plot(x, what = "pairwise", ord = x$working.sequences[[1]]$order$mds$info$locimap$locus)
plot(x, what = "pairwise", ord = rownames(x$working.sequences[[1]]$order$genome$info))
plot(x, what = "pairwise", ord = x$working.sequences[[3]]$order$mds$info$locimap$locus)
plot(x, what = "pairwise", ord = x$working.sequences[[4]]$order$mds$info$locimap$locus)
plot(x, what = "pairwise", ord = x$working.sequences[[5]]$order$mds$info$locimap$locus)
x <- genome_order_per_group(x, gr = 1)
x <- genome_order_per_group(x, gr = 2)
x <- genome_order_per_group(x, gr = 3)
x <- genome_order_per_group(x, gr = 4)
x <- genome_order_per_group(x, gr = 5)
plot(x, what = "pairwise", ord = rownames(x$working.sequences[[1]]$order$genome$info))
plot(x, what = "pairwise", ord = rownames(x$working.sequences[[2]]$order$genome$info))
plot(x, what = "pairwise", ord = rownames(x$working.sequences[[3]]$order$genome$info))
plot(x, what = "pairwise", ord = rownames(x$working.sequences[[4]]$order$genome$info))
plot(x, what = "pairwise", ord = rownames(x$working.sequences[[5]]$order$genome$info))


#### IMPLEMENT SEQUENTIAL ####
# FIXME
### Phasing #####
x <- pairwise_phasing_per_group(x, gr = 1, type = "mds")
x <- pairwise_phasing_per_group(x, gr = 2, type = "mds")
x <- pairwise_phasing_per_group(x, gr = 3, type = "mds")
x <- pairwise_phasing_per_group(x, gr = 4, type = "mds")
x <- pairwise_phasing_per_group(x, gr = 5, type = "mds")
x <- pairwise_phasing_per_group(x, gr = 1, type = "genome")
x <- pairwise_phasing_per_group(x, gr = 2, type = "genome")
x <- pairwise_phasing_per_group(x, gr = 3, type = "genome")
x <- pairwise_phasing_per_group(x, gr = 4, type = "genome")
x <- pairwise_phasing_per_group(x, gr = 5, type = "genome")

# Mapping #
x <- mapping_per_group(x, gr = 1, error = 0.05)
x <- mapping_per_group(x, gr = 2, error = 0.05)
x <- mapping_per_group(x, gr = 3, error = 0.05)
plot_map(x, gr = 1, mrk.names = TRUE)
plot_map(x, gr = 2, mrk.names = TRUE, left.lim = 10, right.lim = 20)
plot_map(x, gr = 3, mrk.names = TRUE)

s.ch1.all <- calc_haplotypes(s.ch1.all)










plot_map(s.ch1.p1, mrk.names = TRUE)
#### Select parent 2####
s.ch1.p2 <- make_sequence(dat, arg = "ch1", info.parent = "p2")
tpt2 <- est_pairwise_rf(s.ch1.p2, ncpus = 7)
# Phasing #
s.ch1.p2 <- pairwise_phasing(input.seq = s.ch1.p2,
                             input.twopt = tpt2,
                             thresh.LOD.ph = 5,
                             max.conf.btnk.p2 = 5)
# Mapping #
s.ch1.p2 <- mapping(input.seq = s.ch1.p2,
                    verbose = TRUE,
                    error = 0.00,# error = 0.0
                    tol = 10e-4)
plot_map(s.ch1.p2, mrk.names = TRUE)
s.ch1.p2 <- mapping(input.seq = s.ch1.p2,
                    verbose = TRUE,
                    error = 0.05,
                    tol = 10e-4)
plot_map(s.ch1.p2, mrk.names = TRUE)
#### Merging maps from p1 and p2 ####
s.ch1.all <- merge_single_parent_maps(input.seq.all = s.ch1.all,
                                      input.seq.p1 = s.ch1.p1,
                                      input.seq.p2 = s.ch1.p2,
                                      input.twopt = tpt.ch1)
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
input.mat <- rf_list_to_matrix(tpt.ch1,
                               thresh.LOD.ph = 5,
                               thresh.LOD.rf = 5,
                               thresh.rf = .5,
                               shared.alleles = TRUE)


plot_map(s.ch1.all, mrk.names = T, xlim = c(-10,105))






#### Two-points ####
tpt <- est_pairwise_rf(s, ncpus = 7)
m <- rf_list_to_matrix(tpt)
plot(m)

#### Grouping ####
lg <- group(m, expected.groups = 5, comp.mat = TRUE, inter = F)

#### CH1 ####
s.ch1.all <- make_sequence(lg, 1, genomic.info = 1)
s.ch1.all
tpt.ch1 <- est_pairwise_rf(s.ch1.all) ## Will remove
#### Select parent 1####
s.ch1.p1 <- make_sequence(s.ch1.all, info.parent = "p1")
s.ch1.p1
# Phasing #
s.ch1.p1 <- pairwise_phasing(input.seq = s.ch1.p1,
                             input.twopt = tpt,
                             thresh.LOD.ph = 3,
                             max.conf.btnk.p1 = 10)
# Mapping #
s.ch1.p1 <- mapping(input.seq = s.ch1.p1,
                    verbose = TRUE,
                    error = 0.00,# error = 0.0
                    tol = 10e-4)
plot_map(s.ch1.p1, mrk.names = TRUE)
#### Select parent 2####
s.ch1.p2 <- make_sequence(s.ch1.all, info.parent = "p2")
s.ch1.p2
# Phasing #
s.ch1.p2 <- pairwise_phasing(input.seq = s.ch1.p2,
                             input.twopt = tpt,
                             thresh.LOD.ph = 3,
                             max.conf.btnk.p1 = 10)
# Mapping #
s.ch1.p2 <- mapping(input.seq = s.ch1.p2,
                    verbose = TRUE,
                    error = 0.00,# error = 0.0
                    tol = 10e-4)
plot_map(s.ch1.p2, mrk.names = TRUE)
#### Merging maps from p1 and p2 ####
s.ch1.all <- merge_single_parent_maps(input.seq.all = s.ch1.all,
                                      input.seq.p1 = s.ch1.p1,
                                      input.seq.p2 = s.ch1.p2,
                                      input.twopt = tpt.ch1)
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
input.mat <- rf_list_to_matrix(tpt.ch1,
                               thresh.LOD.ph = 5,
                               thresh.LOD.rf = 5,
                               thresh.rf = .5,
                               shared.alleles = TRUE)


plot_map(s.ch1.all, mrk.names = T, xlim = c(-10,105))


### Begin
mrk <- "Ch_1_M_1"
drop.seq <- drop_marker(input.seq = s.ch1.all, mrk, reestimate.map = T)
plot_map(drop.seq, 0, 10, mrk.names = T, xlim = c(-2, 11), cex = 1.5)
drop.seq<-calc_haplotypes(drop.seq)
drop.seq <- add_marker(drop.seq, mrk, thresh.rf.to.add = .1)
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


rm(list = ls())
require(mappoly2)
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
                  n.chrom = 12,
                  random = FALSE,
                  seed = 43598)
dat <- read_geno_csv(file.in = "misc/fake_triploid.csv",
                     ploidy.p1 = ploidy.p1,
                     ploidy.p2 = ploidy.p2,
                     name.p1 = "parent_1",
                     name.p2 = "parent_2")
plot(dat)
print(dat, detailed = TRUE)
###Subset test ####
x <- subset(dat)
plot(x)
x1 <- subset(B2721, type = "marker", n = 700)
x1 <- subset(x1, type = "individual", n = 80)
x2 <- subset(B2721, type = "marker", n = 500)
x2 <- subset(x2, type = "individual", n = 100)
x3 <- subset(B2721, type = "marker", n = 400)
x3 <- subset(x3, type = "individual", n = 100)
x <- merge_datasets(x1, x2, x3)
x
plot(x)

####Filter test####
dat <- filter_missing(dat, type = "marker", filter.thres = 0.15, inter = FALSE)
dat <- filter_missing(dat, type = "individual", filter.thres = 0.11, inter = FALSE)
dat <- filter_individuals(dat)
s <- filter_segregation(s, inter = FALSE)
####Two point####
s <- make_sequence(dat, "all")
print(s, detailed = TRUE)
plot(s)
s <- pairwise_rf(s, ncpus = 1)
s
plot(s, type = "rf", fact = 2)
#### Grouping ####
s <- group(s, expected.groups = 12, comp.mat = TRUE, inter = F)
s
plot(s)
#### Ordering ####





#### CH1 ####
s.ch1.p1 <- pairwise_phasing(input.seq = s,
                             lg = 1,
                             info.parent = c("both", "p1", "p2"),
                             thresh.LOD.ph = 3,
                             max.conf.btnk.p1 = 10)




#### Select parent 1####
s.ch1.p1 <- make_sequence(s.ch1.filt, info.parent = "p1")
s.ch1.p1
# Phasing #
s.ch1.p1 <- pairwise_phasing(input.seq = s.ch1.p1,
                             input.twopt = tpt,
                             thresh.LOD.ph = 3,
                             max.conf.btnk.p1 = 10)
# Mapping #
s.ch1.p1 <- mapping(input.seq = s.ch1.p1,
                    verbose = TRUE,
                    error = 0.00,# error = 0.0
                    tol = 10e-4)
plot_map(s.ch1.p1, mrk.names = TRUE)
#### Select parent 2####
s.ch1.p2 <- make_sequence(s.ch1.all, info.parent = "p2")
s.ch1.p2
# Phasing #
s.ch1.p2 <- pairwise_phasing(input.seq = s.ch1.p2,
                             input.twopt = tpt,
                             thresh.LOD.ph = 3,
                             max.conf.btnk.p1 = 10)
# Mapping #
s.ch1.p2 <- mapping(input.seq = s.ch1.p2,
                    verbose = TRUE,
                    error = 0.00,# error = 0.0
                    tol = 10e-4)
plot_map(s.ch1.p2, mrk.names = TRUE)
#### Merging maps from p1 and p2 ####
s.ch1.all <- merge_single_parent_maps(input.seq.all = s.ch1.all,
                                      input.seq.p1 = s.ch1.p1,
                                      input.seq.p2 = s.ch1.p2,
                                      input.twopt = tpt.ch1)
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
input.mat <- rf_list_to_matrix(tpt.ch1,
                               thresh.LOD.ph = 5,
                               thresh.LOD.rf = 5,
                               thresh.rf = .5,
                               shared.alleles = TRUE)
plot_map(s.ch1.all, mrk.names = T, xlim = c(-10,105))



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





rm(list = ls())
require(mappoly2)
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
                  n.chrom = 12,
                  random = FALSE,
                  seed = 43598)
dat <- read_geno_csv(file.in = "misc/fake_triploid.csv",
                     ploidy.p1 = ploidy.p1,
                     ploidy.p2 = ploidy.p2,
                     name.p1 = "parent_1",
                     name.p2 = "parent_2")
plot(dat)
print(dat, detailed = TRUE)
###Subset test ####
x <- subset(dat)
plot(x)
x1 <- subset(B2721, type = "marker", n = 700)
x1 <- subset(x1, type = "individual", n = 80)
x2 <- subset(B2721, type = "marker", n = 500)
x2 <- subset(x2, type = "individual", n = 100)
x3 <- subset(B2721, type = "marker", n = 400)
x3 <- subset(x3, type = "individual", n = 100)
x <- merge_datasets(x1, x2, x3)
x
plot(x)

####Filter test####
dat <- filter_missing(dat, type = "marker", filter.thres = 0.15, inter = FALSE)
dat <- filter_missing(dat, type = "individual", filter.thres = 0.11, inter = FALSE)
dat <- filter_individuals(dat)
s <- filter_segregation(s, inter = FALSE)
####Two point####
s <- make_sequence(dat, "all")
print(s, detailed = TRUE)
plot(s)
s <- pairwise_rf(s, ncpus = 1)

B2721 <- read_geno_csv(file.in = "inst/extdata/B2721.csv",
                       ploidy.p1 = 4,
                       ploidy.p2 = 4,
                       name.p1 = B2721$name.p1,
                       name.p2 = B2721$name.p2)
save(B2721, file = "data/example_data.rda")
