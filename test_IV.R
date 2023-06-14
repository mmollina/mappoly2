###### OLD ####
rm(list = ls())
require(mappoly2)
source("~/repos/official_repos/misc/simulation.R")
ploidy.p1 = 2
ploidy.p2 = 4
n.mrk <- 300
ph<-test_simulate(ploidy.p1 = ploidy.p1,
                  ploidy.p2 = ploidy.p2,
                  fpath = "~/repos/official_repos/misc/fake_triploid.csv",
                  n.mrk = n.mrk,
                  n.ind = 200,
                  map.length =100,
                  miss.perc = 10,
                  n.chrom = 1,
                  random = FALSE,
                  seed = 2986876)
####Read ####
dat <- read_geno_csv(file.in = "~/repos/official_repos/misc/fake_triploid.csv",
                     ploidy.p1 = ploidy.p1, ploidy.p2 = ploidy.p2)
dat
#### Select parent ####
s <- make_sequence(dat, "all", info.parent = "both")
#s <- make_sequence(dat, "all", info.parent = "p1")
#s <- make_sequence(dat, "all", info.parent = "p2")
#### Two-points ####
tpt <- est_pairwise_rf(s, ncpus = 7)

#### Phasing ####
s <- pairwise_phasing(input.seq = s,
                      input.twopt = tpt,
                      thresh.LOD.ph = 15,
                      max.conf.btnk.p1 = 1)
s
#### Mapping ####
s1 <- mapping(input.seq = s,
              verbose = TRUE,
              error = 0.00,# error = 0.0
              tol = 10e-4)

plot(cumsum(imf_h(s1$phases[[1]]$rf)))

x<- NULL
for(i in seq(30,0,-5)){
  m <- rf_list_to_matrix(tpt,
                         thresh.LOD.ph = i,
                         shared.alleles = TRUE)
  mrk.pos <- rownames(s1$phases[[1]]$p2)
  mrk.id <- setdiff(s1$mrk.names, mrk.pos)
  InitPh <- s1$phases[[1]]$p2
  dose.vec <- s1$data$dosage.p2[mrk.id]
  S <- m$Sh.p2[mrk.id, mrk.pos]
  verbose <- TRUE
  L <- mappoly2:::phasing_one(mrk.id, dose.vec, S, InitPh, verbose)

  L

  x <- cbind(x, sort(sapply(L, nrow)))
}
image(x)

########## Testing vs_inserted_mrk #########
i<-5
m <- rf_list_to_matrix(tpt,
                       thresh.LOD.ph = i,
                       shared.alleles = TRUE)
mrk.pos <- rownames(s1$phases[[1]]$p2)
mrk.id <- setdiff(s1$mrk.names, mrk.pos)
##P1
InitPh <- s1$phases[[1]]$p1
dose.vec <- s1$data$dosage.p1[mrk.id]
S <- m$Sh.p1[mrk.id, mrk.pos]
verbose <- TRUE
L1 <- mappoly2:::phasing_one(mrk.id, dose.vec, S, InitPh, verbose)
L1[[48]]
##P2
InitPh <- s1$phases[[1]]$p2
dose.vec <- s1$data$dosage.p2[mrk.id]
S <- m$Sh.p2[mrk.id, mrk.pos]
verbose <- TRUE
L2 <- mappoly2:::phasing_one(mrk.id, dose.vec, S, InitPh, verbose)
L2[[48]]
L1[[48]]
G <- s1$data$geno.dose[mrk.id[48], ,drop = TRUE]
pedigree <- matrix(rep(c(1,
                         2,
                         s1$data$ploidy.p1,
                         s1$data$ploidy.p2, 1),
                       s1$data$n.ind),
                   nrow = s1$data$n.ind,
                   byrow = TRUE)
s1<-calc_haplotypes(s1)
mrk.id[48]
s1$phases[[1]]$p1[18:19,]
s1$phases[[1]]$p2[18:19,]
homolog_prob <- as.matrix(s1$phases[[1]]$haploprob[,18:19+2])

L2[[48]]
L1[[48]]
PH1 <- list(L1[[48]][1,,drop=TRUE], L2[[48]][1,,drop = TRUE])
z1<-mappoly2:::est_hmm_map_biallelic_insert_marker(PH1,
                                                   G,
                                                   pedigree,
                                                   homolog_prob,
                                                  rf = c(0.01,0.01),
                                                  verbose = TRUE,
                                                  detailed_verbose = FALSE,
                                                  tol = 10e-4,
                                                  ret_H0 = FALSE)
PH2 <- list(L1[[48]][1,,drop=TRUE], L2[[48]][2,,drop = TRUE])
z2<-mappoly2:::est_hmm_map_biallelic_insert_marker(PH2,
                                                   G,
                                                   pedigree,
                                                   homolog_prob,
                                                   rf = c(0.01,0.01),
                                                   verbose = TRUE,
                                                   detailed_verbose = FALSE,
                                                   tol = 10e-4,
                                                   ret_H0 = FALSE)
z1
z2

