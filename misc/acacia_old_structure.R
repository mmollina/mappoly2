rm(list = ls())
rm(list = ls())
require(mappoly2)
P1 <- "AC46"
P2 <- "AC78"
dat <- read_geno_csv(file.in = paste0("~/repos/collaborations/acacia/full_sib_data/mappoly2_", P1, "-", P2,".csv"),
                   ploidy.p1 = 2,
                   ploidy.p2 = 2,
                   name.p1 = P1,
                   name.p2 = P2)
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
s <- make_sequence(dat, "ch1")
s <- filter_segregation(s)

####Two point####
print(s, detailed = TRUE)
plot(s)
s <- pairwise_rf(s, ncpus = 8)
s
plot(s, type = "rf", fact = 2)
#### Grouping ####
s <- group(s, expected.groups = 1, comp.mat = TRUE, inter = F)
s
plot(s)
#### Ordering ####



#### Select parent 1####
s.ch1.p1 <- make_sequence(s, info.parent = "p1")
s.ch1.p1
# Phasing #
s.ch1.p1 <- pairwise_phasing(input.seq = s.ch1.p1,
                             input.twopt = s$pairwise,
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





