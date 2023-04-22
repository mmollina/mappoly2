rm(list = ls())
require(mappoly2)
file.in <- "~/repos/official_repos/mappoly2/inst/extdata/B2721.csv"
input.data <- read_geno_csv(file.in  = file.in,
                            ploidy.p1 = 4,
                            ploidy.p2 = 4,
                            name.p1 = "Atlantic",
                            name.p2 = "B1829-5", filter.non.conforming = T, filter.redundant = F)
#save(B2721, file = "./data/example_data.rda")
input.data
plot(input.data)
input.seq <- y1 <- make_seq(input.data, "all")
y2 <- make_seq(input.data, c("ch1", "CHR4"))
y3 <- make_seq(input.data, arg = y1$mrk.names[1:100])
y4 <- make_seq(input.data, 1:100)
tpt1 <- est_pairwise_rf(y4, ncpus = 2)
tpt2 <- est_pairwise_rf(y4, ncpus = 1)
identical(tpt1, tpt2)
identical(y3,y4)
y1;y2;y3;y4

dat.filt <- filter_missing(input.data)
dat.filt <- filter_missing(dat.filt)
dat.filt <- filter_missing(dat.filt, type = "individual", filter.thres = 0.04)
dat.filt <- filter_individuals(dat.filt)
dat.filt <- filter_individuals(dat.filt)
sf1 <- filter_segregation(dat.filt)
plot(sf1)
s2 <- make_seq(dat.filt, "all")
sf2 <- filter_segregation(input.obj = s2)
y5 <- make_seq(sf1, "ch1")
plot(y5)
mrk.id<-names(sort(y5$data$genome.pos[y5$mrk.names]))
tpt5 <- est_pairwise_rf(y5, ncpus = 7)
a<-rf_list_to_matrix(input.twopt = tpt5)
plot(a, ord = mrk.id)


source("~/repos/official_repos/mappoly2/misc/simulation.R")
for(ploidy.p1 in c(2,4,6)){
  for(ploidy.p2 in c(2,4,6)){
    ph<-test_simulate(ploidy.p1 = ploidy.p1,
                  ploidy.p2 = ploidy.p2,
                  fpath = "./misc/test_sim_data.csv",
                  n.mrk = 100,
                  n.ind = 3000,
                  map.length = 100)
    dat <- read_geno_csv(file.in = "./misc/test_sim_data.csv", ploidy.p1 = ploidy.p1, ploidy.p2 = ploidy.p2)
    dat
    input.seq<-make_seq(dat, "all")
    tpt1<-est_pairwise_rf(input.seq, ncpus = 1)
    m<-rf_list_to_matrix(tpt1)
    plot(m, main.text = paste(ploidy.p1, "--x--", ploidy.p2, "---> swap: ", dat$swap.parents), ord = 35:65, index = T)
    readline("next...")
  }
}
