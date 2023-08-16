rm(list = ls())
require(mappoly2)
#### Read/QAQC ####
x <- read_geno_csv(file.in = "~/repos/collaborations/acacia/full_sib_data/mappoly2_AC48-AC78.csv",
                     ploidy.p1 = 2,
                     ploidy.p2 = 2,
                     name.p1 = "AC48",
                     name.p2 = "AC78")
plot(x)
class(x)
x
x <- filter_data(x, mrk.thresh = .08, ind.thresh = .7)
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
x <- set_initial_sequence(x, c("ch1", "ch4", "ch6"))
plot(x, what = "initiated")
#### Compute RF ####
x <- pairwise_rf(x, ncpus = 8)
plot(x, what = "pairwise", fact = 10)
plot(x, what = "pairwise", type = "lod", fact = 10)
### FIXME: rf_filter for initial vs. rf_filter for working
### Make two functions? rf_filter_initial and rf_filter_working ????
x <- rf_filter(x, 5, 5, .1, c(0.01, 1))

#### Grouping ####
x <- group(x, expected.groups = 3)
plot(x, what = "group")
x
#### Set Working Sequence ####
mrks.ch4 <- mappoly2:::get_markers_from_grouped_and_chromosome(x, lg = 2, ch = "ch4")






x <- rf_filter(x, 5, 5, .1, c(0.01, 1))


plot(x, what = "rf_filtered")
x <- group(x, expected.groups = 13, comp.mat = T)
x
plot_mappoly2_group(x$linkage.groups)

## check print for each step
## check plot for each step
## implement make_sequence






x <- group(x, expected.groups = 5, comp.mat = TRUE, inter = F)









#### Select parent 1####
s.ch1.p1 <- make_sequence(dat, arg = "ch1", info.parent = "p1")
s.ch1.p1
tpt1 <- est_pairwise_rf(s.ch1.p1, ncpus = 7)
# Phasing #
s.ch1.p1 <- pairwise_phasing(input.seq = s.ch1.p1,
                             input.twopt = tpt1,
                             thresh.LOD.ph = 5,
                             max.conf.btnk.p1 = 10)
# Mapping #
s.ch1.p1 <- mapping(input.seq = s.ch1.p1,
                    verbose = TRUE,
                    error = 0.00,# error = 0.0
                    tol = 10e-4)
plot_map(s.ch1.p1, mrk.names = TRUE)
s.ch1.p1 <- mapping(input.seq = s.ch1.p1,
                    verbose = TRUE,
                    error = 0.05,
                    tol = 10e-4)
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





