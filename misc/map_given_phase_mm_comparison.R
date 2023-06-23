rm(list = ls())
require(mappoly2)
require(crayon)
source("~/repos/official_repos/misc/simulation.R")
ploidy.p1 = 4
ploidy.p2 = 4
n.mrk <- 100
plot(c(1,n.mrk), c(0, 3), type = "n")
time_mp2 <- numeric(100)
time_mp <- numeric(100)
abline(h = 0, lty = 2)
for(i in 1:100){
  ph<-test_simulate(ploidy.p1 = ploidy.p1,
                    ploidy.p2 = ploidy.p2,
                    fpath = "~/repos/official_repos/misc/fake_triploid.csv",
                    n.mrk = n.mrk,
                    n.ind = 100,
                    map.length = 150,
                    miss.perc = 0,
                    n.chrom = 1,
                    random = FALSE, seed = i)
  dat <- read_geno_csv(file.in = "~/repos/official_repos/misc/fake_triploid.csv", ploidy.p1 = ploidy.p1, ploidy.p2 = ploidy.p2, filter.non.conforming = F, filter.redundant = F)
  s1 <- make_sequence(dat, "all")
  id <- s1$mrk.names
  PH <- list(ph[[1]]$p1[id,], ph[[1]]$p2[id,])
    time_mp2[i] <- system.time(map <- hmm_map_reconstruction(s1,
                                                           PH,
                                                           verbose = TRUE,
                                                           tol = 10e-4))[3]
  #####MAPPOLY1####
  dat1 <- read.csv("~/repos/official_repos/misc/fake_triploid.csv")
  dat1 <- mappoly::table_to_mappoly(dat1[, -c(6,7)], ploidy = 4, filter.non.conforming = F, elim.redundant = F)
  smp <- mappoly::make_seq_mappoly(dat1, 1:300)
  PH2 <- list(P = mappoly::ph_matrix_to_list(PH[[1]]),
              Q = mappoly::ph_matrix_to_list(PH[[2]]))
  time_mp[i] <- system.time(maps1 <- mappoly::est_rf_hmm_single(input.seq = smp,
                                                                input.ph.single = PH2,
                                                                tol = 10e-4,
                                                                high.prec = FALSE,
                                                                verbose = T))[3]
  lines(rep(i,2), c(time_mp[i], time_mp2[i]))
  points(rep(i,2), c(time_mp[i], time_mp2[i]),  col = c(2, 4), pch = 19)
  points(i, map$loglike/maps1$loglike, col = 3, pch = 3)
  points(i, sum(maps1$seq.rf -  map$rec.frac))
}

