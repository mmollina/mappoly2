rm(list = ls())
require(mappoly2)
source("~/repos/official_repos/misc/simulation.R")
ploidy.p1 = 2
ploidy.p2 = 4
ph<-test_simulate(ploidy.p1 = ploidy.p1,
                  ploidy.p2 = ploidy.p2,
                  fpath = "~/repos/official_repos/misc/fake_triploid.csv",
                  n.mrk = 100,
                  n.ind = 100,
                  map.length = 10,
                  miss.perc = 0,
                  n.chrom = 3,
                  random = FALSE,
                  seed = 1234)
dat <- read_geno_csv(file.in = "~/repos/official_repos/misc/fake_triploid.csv", ploidy.p1 = ploidy.p1, ploidy.p2 = ploidy.p2, filter.redundant = F)
s <- make_sequence(dat, )
tpt <- est_pairwise_rf(s, ncpus = 8)
m <- rf_list_to_matrix(tpt)
gen.ord <- get_genome_order(input.seq = s)
plot(m)
plot(m, ord = rownames(gen.ord$ord))
#s1 <- make_sequence(dat, 1:2)
#tpt1<-est_pairwise_rf(s1)
#ph[[1]]$p1[1:4,]
#ph[[1]]$p2[1:4,]
#tpt1$pairwise$`1-2`
#dat$chrom[sample(1:900, 100)] <- NA
plot(dat)
s1 <- make_sequence(dat, "all")
s1
plot(s1)
s.ch1 <- make_sequence(dat, "ch1")
s.ch1
plot(s.ch1)
s.ch2.3 <- make_sequence(dat, c("ch2", "ch3"))
plot(s.ch2.3)
s.ch2 <- make_sequence(input.obj = s.ch2.3, arg = "ch2")
s.ch2
dat.filt <- filter_missing(dat, filter.thres = 0.08, inter = F)
dat.filt <- filter_missing(dat.filt)
dat.filt <- filter_missing(dat.filt, type = "individual", filter.thres = 0.07)
dat.filt <- filter_individuals(dat.filt)
dat.filt <- filter_individuals(dat.filt)
sf1 <- filter_segregation(input.obj = dat.filt, inter = FALSE)
plot(sf1)
s.ch1.2.3.filt<- make_sequence(dat.filt, "all")
tpt <- est_pairwise_rf(s.ch1.2.3.filt, ncpus = 4)
m <- rf_list_to_matrix(tpt)
gen.ord <- get_genome_order(input.seq = s.ch1.2.3.filt)
plot(m)
plot(m, ord = rownames(gen.ord$ord))

s.ch1.2.3.rf_filt<-rf_snp_filter(tpt, probs = c(0.01,1))
plot(m)
plot(m, ord = c(20:30), index = T)
gr <- group(m,expected.groups = 3,comp.mat = TRUE, inter = FALSE)
gr
plot(gr)
s1 <- make_sequence(input.obj = gr, arg = 1)
s1 <- make_sequence(dat, arg = "ch1")
tpt1 <- est_pairwise_rf(s1, ncpus = 8)
m1 <- rf_list_to_matrix(tpt1)
o1 <- mds(m1)
so1 <- make_sequence(o1)
plot(so1$data$genome.pos[so1$mrk.names])
plot(m1, ord = so1)

## constructing class phased sequence





rm(list = ls())
require(mappoly2)
source("~/repos/official_repos/misc/simulation.R")
ploidy.p1 = 2
ploidy.p2 = 4
ph<-test_simulate(ploidy.p1 = ploidy.p1,
                  ploidy.p2 = ploidy.p2,
                  fpath = "~/repos/official_repos/misc/fake_triploid.csv",
                  n.mrk = 300,
                  n.ind = 100,
                  map.length = 100,
                  miss.perc = 0,
                  n.chrom = 3,
                  random = FALSE,
                  seed = 1234)
dat <- read_geno_csv(file.in = "~/repos/official_repos/misc/fake_triploid.csv", ploidy.p1 = ploidy.p1, ploidy.p2 = ploidy.p2, filter.redundant = F)
s1<-make_sequence(dat, 1:20)
id <- s1$mrk.names
G <- s1$data$geno.dose[id, ]
PH <- list(ph[[1]]$p1[id,], ph[[1]]$p2[id,])
pedigree <- matrix(rep(c(1,
                         2,
                         s1$data$ploidy.p1,
                         s1$data$ploidy.p2, 1),
                       s1$data$n.ind),
                   nrow = s1$data$n.ind,
                   byrow = TRUE)
head(pedigree)

bla<-mappoly2:::est_hmm_map_biallelic(PH = PH,
                                 G = G,
                                 pedigree = pedigree,
                                 verbose = TRUE,
                                 tol = 10e-4,
                                 ret_H0 = FALSE)


