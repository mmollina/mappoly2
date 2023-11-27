rm(list = ls())
require(mappoly2)
require(tidyverse)
dat <- read_geno_csv(file.in = "~/repos/BI_sweetpotato.csv",
                     ploidy.p1 = 6,
                     ploidy.p2 = 6,
                     name.p1 = "Beauregard",
                     name.p2 = "Regal")
dat
plot(dat)
dat <- filter_missing(dat, type = "marker", filter.thres = 0.10, inter = F)
dat <- filter_missing(dat, type = "individual", filter.thres = 0.10, inter = F)
dat <- filter_individuals(dat)


s <- filter_segregation(dat, inter = FALSE)
print(s, detailed = TRUE)
plot(s)

#### Two-points ####
tpt <- est_pairwise_rf(s, ncpus = 7)
m <- rf_list_to_matrix(tpt, thresh.LOD.ph = 1, thresh.LOD.rf = 1)
plot(m)

so <- get_genome_order(s)
plot(so)
so <- make_sequence(so)
plot(m, ord = so, fact = 5)

#### Grouping ####
gr <- group(m, expected.groups = 20, inter = FALSE, comp.mat = TRUE)
plot(gr)
gr
heatmap(gr$seq.vs.grouped.snp, Colv = NA, xlab = "Chrom",  ylab = "LG")
## convert to tibble, add row identifier, and shape "long"
dat2 <-
  gr$seq.vs.grouped.snp %>%
  as_tibble() %>%
  rownames_to_column("Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value") %>%
  rename(Linkage_group = Var1) %>% rename(Chrom = Var2)

ggplot(dat2, aes(Linkage_group, Chrom)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = round(value, 1))) +
  scale_fill_gradient(low = "white", high = "red")

#### CH1 ####
dat.sim <- read_geno_csv(file.in = "misc/fake_triploid.csv",
                         ploidy.p1 = ploidy.p1,
                         ploidy.p2 = ploidy.p2,
                         name.p1 = "parent_1",
                         name.p2 = "parent_2")
s <- make_sequence(dat.sim, arg = "ch1")
tpt <- est_pairwise_rf(s) ## Will remove
y1<-numeric(100)
x1<- seq(0, 20, length.out = 100)
for(i in 1:100){
  m <- rf_list_to_matrix(tpt, thresh.LOD.ph =x1[i])
  y1[i] <- sum(!is.na(m$lod.mat))/length(m$lod.mat)
}
#plot(y1 ~ x1, xlab = "LOD", ylab = "percetange filled", xlim = c(-1, max(x1)), ylim = c(0,1), type = "l")
points(y1 ~ x1, col = 7, type = "l", lwd = 2)





s.ch1.all <- make_sequence(dat, arg = "ch1", genomic.info = 1)
s.ch1.all
tpt.ch1 <- est_pairwise_rf(s.ch1.all) ## Will remove
y2<-numeric(100)
x2<- seq(0, 20, length.out = 100)
for(i in 1:100){
  m.ch1 <- rf_list_to_matrix(tpt.ch1, thresh.LOD.ph = x2[i])
  y2[i] <- sum(!is.na(m.ch1$lod.mat))/length(m.ch1$lod.mat)
}
points(y2 ~ x2, col = 2, type = "l")
for(i in 1:length(x1))
lines(c(x1[i],x2[i]), c(y1[i], y2[i]), col = 3)







s.ch1 <- pairwise_phasing(input.seq = s.ch1.all,
                          input.twopt = tpt.ch1,
                          thresh.LOD.ph = 1,
                          max.conf.btnk.p1 = 20)



#### Select parent 1####
s.ch1.p1 <- make_sequence(s.ch1.all, info.parent = "p1")
s.ch1.p1
# Phasing #
s.ch1.p1 <- pairwise_phasing(input.seq = s.ch1.p1,
                             input.twopt = tpt.ch1,
                             thresh.LOD.ph = .5,
                             max.conf.btnk.p1 = 200)
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
                             thresh.LOD.ph = .5,
                             max.conf.btnk.p1 = 200)
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
plot_map(s.ch1.all, 60, 80, mrk.names = T, xlim = c(58, 81), cex = 1.5)
x11()
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





