rm(list = ls())
require(mappoly2)
require(crayon)
source("~/repos/official_repos/misc/twopt_phasing.R")
source("~/repos/official_repos/misc/simulation.R")
ph<-test_simulate(ploidy.p1 = 6,
                  ploidy.p2 = 4,
                  fpath = "~/repos/official_repos/misc/fake_triploid.csv",
                  n.mrk = 1000,
                  n.ind = 1000,
                  map.length = 100,
                  miss.perc = 0,
                  n.chrom = 1,
                  random = FALSE,
                  seed = 3489193)
dat <- read_geno_csv(file.in = "~/repos/official_repos/misc/fake_triploid.csv",
                     ploidy.p1 = ncol(ph[[1]]$p1), ploidy.p2 = ncol(ph[[1]]$p2))
s1 <- make_sequence(dat, "all")
mrk.id <- s1$mrk.names
PH <- list(ph[[1]]$p1[mrk.id,], ph[[1]]$p2[mrk.id,])
tpt <- est_pairwise_rf(s1, ncpus = 6)
w <- sapply(tpt$pairwise, rownames)
x <- t(sapply(w, function(y) apply(stringr::str_split_fixed(y, "-", 2), 2, function(x) max(as.numeric(x)))))
apply(x, 2, range)
m <- rf_list_to_matrix(tpt, shared.alleles = T, thresh.LOD.ph = 1)
plot(m)
S <- m$Sh.p1[mrk.id, mrk.id]
dose <- s1$data$dosage.p1[mrk.id]
ploidy <- s1$data$ploidy.p1
t2 <- system.time(H2 <- mappoly2:::twopt_phasing_cpp(mrk.id,ploidy,dose,S,max_conf_number = 10))



length(mrk.id)
dose <- s1$data$dosage.p2[mrk.id]
ploidy <- s1$data$ploidy.p2
S <- m$Sh.p2[mrk.id, mrk.id]
range(m$Sh.p2[mrk.id, mrk.id], na.rm = T)
range(m$Sh.p1[mrk.id, mrk.id], na.rm = T)
image(S)
dim(S)
Ph <- PH[[1]][mrk.id,]

{cat(bgRed('~~~~~~~~~~~~~~~~~~R~~~~~~~~~~~~~~~~~~\n'))
  # #################R implementation#########################
  t1 <- system.time(H1 <- twopt_phasing(mrk.id, dose, ploidy, S, Ph))
  length(H1)
  Ph <- PH[[1]][mrk.id,]
  for(i in 1:length(H1)){
    p1 <- mappoly::ph_matrix_to_list(H1[[i]])
    p2 <- mappoly::ph_matrix_to_list(Ph[rownames(H1[[i]]),])
    a<-mappoly::compare_haplotypes(dat$ploidy.p1, p1, p2)
    if(a$is.same.haplo){
      cat(bgGreen('\nCorrect!\n'))
      cat(bgCyan('haplo ord?: ', paste0(a$haplo.ord, collapse = "-"), '\n'))
    }
  }
  cat(bgYellow('elapsed time:',round(t1[3],4),'\n'))
  cat(bgYellow('N Config.',length(H1),'\n'))
  for(i in 1:length(H1)){
    A <- tcrossprod(H1[[i]])
    res <- S[colnames(A), colnames(A)] - A
    x1 <- sum(res, na.rm = T)
    cat(bgGreen('Check (0=OK) --> ',  x1,'\n'))
  }
  cat(bgBlue('N. mrk --> ',  nrow(H1[[1]]),'\n'))
  cat(bgRed('~~~~~~~~~~~~~~~~~~Rcpp~~~~~~~~~~~~~~~~~~\n'))
  #################RCpp implementation#########################

  n.mrk <- n.conf <- NULL
x <- seq(0,5,2)
  for(j in x){
    m <- rf_list_to_matrix(tpt, shared.alleles = T, thresh.LOD.ph = j)
    S <- m$Sh.p1[mrk.id, mrk.id]
    t2 <- system.time(H2 <- mappoly2:::twopt_phasing_cpp(mrk.id,ploidy,dose,S,max_conf_number = 300000))
    length(H2$phase_configs)
    for(i in 1:length(H2$phase_configs))
      rownames(H2$phase_configs[[i]]) <- H2$marker_names
    for(i in 1:length(H2$phase_configs))
    {
      p1 <- mappoly::ph_matrix_to_list(H2$phase_configs[[i]])
      p2 <- mappoly::ph_matrix_to_list(Ph[rownames(H2$phase_configs[[i]]),])
      a<-mappoly::compare_haplotypes(dat$ploidy.p1, p1, p2)
      if(a$is.same.haplo){
        cat(bgGreen('\nCorrect!\n'))
        cat(bgCyan('haplo ord?: ', paste0(a$haplo.ord, collapse = "-"), '\n'))
      }
    }
    n.conf <- c(n.conf, length(H2$phase_configs))
    cat(bgYellow('elapsed time:',round(t2[3],4),'\n'))
    cat(bgYellow('N Config.',length(H2$phase_configs),'\n'))
    for(i in 1:length(H2$phase_configs)){
      A <- tcrossprod(H2$phase_configs[[i]])
      res <- S[colnames(A), colnames(A)] - A
      x2 <- sum(res, na.rm = T)
      cat(bgGreen('Check (0=OK) --> ',  x2,'\n'))
    }
    cat(bgBlue('N. mrk --> ',  length(H2$marker_names),'\n'))
    n.mrk <- c(n.mrk, length(H2$marker_names))
  }
}
