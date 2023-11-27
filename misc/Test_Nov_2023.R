rm(list = ls())
require(mappoly2)
setwd("~/repos/official_repos/mappoly2/")
source("misc/simulation.R")
ploidy.p1 = 4
ploidy.p2 = 4
n.mrk <- 300
map.length = 100
ph<-test_simulate(ploidy.p1 = ploidy.p1,
                  ploidy.p2 = ploidy.p2,
                  fpath = "misc/fake_data.csv",
                  n.mrk = n.mrk,
                  n.ind = 300,
                  map.length = map.length,
                  miss.perc = 5,
                  n.chrom = 3,
                  random = FALSE,
                  seed = 2986)
####Reading Data####
x <- read_geno_csv(file.in = "misc/fake_data.csv",
                   ploidy.p1 = ploidy.p1,
                   ploidy.p2 = ploidy.p2)









#### Initial QA/QC ####
x <- filter_data(x, mrk.thresh = .08, ind.thresh = .7)
plot(x)
plot(x, 'screened')
x <- filter_individuals(x)
plot(x, 'screened')

#### Set Initial Sequence ####
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





