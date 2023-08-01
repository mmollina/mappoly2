ploidy  = 6
n.chr = 1
#### Functions ####
require(tidyverse)
sim_map  <- function(ploidy = 4,                     # ploidy level
                     map.length = 100,                 # distance in cM
                     n.ind = 200 ,                     # number of individuals
                     n.mrk = 200,                      # number of markers
                     prob.dose = c(0.45,0.45,0.1),
                     seed.for.config = 40,             # seed for simulating the linkage phase configuration
                     seed.for.pop = 1,
                     prefPairing  =  0,                # preferential pairing
                     quadrivalents  =  0,              # quadrivalent formation rate
                     centromere  =  20,                # long arm 4 times as long as the short arm
                     phase.number.limit = +Inf,        # Limit of number of phases to test in each marker addition 
                     sub.map.size.diff.limit = +Inf,   # Limit of map increment in each marker addition
                     natural.pairing  =  0)            #probabilities of bivalents and quadrivalents are determined by the “quadrivalents” argument
{
  genome.pos <- chrom <- w <- NULL
  for(i in 1:n.chr)
  {
    ph.temp<-mappoly::sim_homologous(ploidy = ploidy,
                                     n.mrk=n.mrk,
                                     max.d=ploidy/2,    #max dosage number
                                     max.ph=ploidy/2,
                                     restriction = FALSE,
                                     prob.dose = prob.dose,
                                     seed = seed.for.config)
    php<-mappoly::ph_list_to_matrix(ph.temp$hom.allele.p, ploidy)
    phq<-mappoly::ph_list_to_matrix(ph.temp$hom.allele.q, ploidy)
    w<-rbind(w, cbind(php,phq))
    chrom <- c(chrom, rep(i, n.mrk))
    genome.pos <- c(genome.pos, 1:n.mrk)
  }
  file.prefix <- paste0("sim_ploidy_", ploidy,
                        "seed_ph", seed.for.config,
                        "seed_pop", seed.for.pop,
                        "_pp_", prefPairing,
                        "_qua_", quadrivalents,
                        "_nch_", n.chr,
                        paste0(prob.dose, collapse = "-"),
                        "number_", i)
  #########PedigreeSim simulation
  #Simulation parameters
  write(x = paste("PLOIDY =", ploidy), file = paste0(file.prefix, ".par"))
  write(x = "MAPFUNCTION = HALDANE", file = paste0(file.prefix, ".par"), append = TRUE)
  write(x = "MISSING = NA", file = paste0(file.prefix, ".par"), append = TRUE)
  write(x = "POPTYPE = F1", file = paste0(file.prefix, ".par"), append = TRUE)
  write(x = paste("POPSIZE =", n.ind), file = paste0(file.prefix, ".par"), append = TRUE)
  write(x = paste("NATURALPAIRING =", natural.pairing), file = paste0(file.prefix, ".par"), append = TRUE)
  write(x = paste0("FOUNDERFILE = ", file.prefix, ".gen"), file = paste0(file.prefix, ".par"), append = TRUE)
  write(x = paste0("CHROMFILE = ", file.prefix, ".chrom"), file = paste0(file.prefix, ".par"), append = TRUE)
  write(x = paste0("MAPFILE = ", file.prefix, ".map"), file = paste0(file.prefix, ".par"), append = TRUE)
  write(x = paste0("OUTPUT = ", file.prefix, "_out"), file = paste0(file.prefix, ".par"), append = TRUE)
  write(x = paste("SEED =", seed.for.pop), file = paste0(file.prefix, ".par"), append = TRUE)
  #chromosome file
  x <- a <- NULL
  for(j in 1:n.chr){
    a <- rbind(a, c(paste0("Ch",j), map.length,	centromere,		prefPairing,		quadrivalents)) 
    x <- rbind(x, data.frame(paste0("M",(1:n.mrk) + n.mrk * (j-1)), paste0("Ch",j), seq(0, map.length,length.out = n.mrk)))
  }
  write.table(a, file =paste0(file.prefix, ".chrom"), quote = FALSE, row.names = FALSE, col.names = c("chromosome",	"length",	"centromere",	"prefPairing", "quadrivalents"))
  write.table(x, file =paste0(file.prefix, ".map"), quote = FALSE, row.names = FALSE, col.names = c("marker", "chromosome", "position"))
  ## parent genotypes
  y<-data.frame(as.character(x[,1]), w)    
  write.table(y, file =paste0(file.prefix, ".gen"), quote = FALSE, row.names = FALSE, col.names = c("marker", paste0("P1_", 1:ploidy), paste0("P2_", 1:ploidy)))
  ## PedigreeSim
  system(paste0("java -jar PedigreeSim.jar ", file.prefix, ".par"), intern=TRUE)
  ############End of simulation
  ## Reading data from pedigreesim
  dat.g<-read.table(file = paste0(file.prefix, "_out_genotypes.dat"), header = TRUE, row.names = 1)
  system(paste0("rm ", file.prefix, "*"), intern=TRUE)
  G<-NULL
  ct<-1
  for(j in 1:(n.ind+2))
  {
    dat.g[,ct:(ct+ploidy-1)]
    G<-cbind(G, apply(dat.g[,ct:(ct+ploidy-1)], 1, function(x) sum(x==1)))
    ct<-ct+ploidy
  }
  P<-G[,c(1:2)]
  colnames(P) <- c("P1", "P2")
  G<-G[,-c(1:2)]
  colnames(G)<-paste0("Ind_", 1:n.ind)
  ref <- sample(c("A", "T", "C", "G"), replace = T, size = nrow(P))
  alt <- character(length(ref))
  for(i in 1:length(ref))
    alt[i] <- sample(setdiff(c("A", "T", "C", "G"), ref[i]), size = 1)
  DF <- data.frame(snp_id = rownames(P),
                   P,
                   chrom = chrom,
                   genome_pos = genome.pos,
                   ref = ref,
                   alt = alt,
                   G)
  dat <- mappoly2::table_to_mappoly(DF, ploidy.p1 = ploidy, ploidy.p2 = ploidy, verbose = FALSE)
  return(dat)
}

#### Sim ####
res <- NULL
require(mappoly2)
setwd("~/repos/Autopolyploid_Linkage/src/simulation2/")
for(n.ind in c(50, 100, 150, 200, 250, 300)){
  cat("ploidy: ", ploidy, " / n.ind: ", n.ind)
  dat.sim <- sim_map(ploidy = ploidy, 
                     map.length = 100, 
                     n.ind = n.ind, 
                     n.mrk = 300, 
                     prob.dose = c(0.4, 0.4, 0.1, 0.1)[1:(1 + ploidy/2)])
  s <- make_sequence(dat.sim, arg = "ch1")
  tpt <- est_pairwise_rf(s, verbose = FALSE) ## Will remove
  y1<-numeric(40)
  x1<- seq(0, 40, length.out = 40)
  for(i in 1:40){
    m <- rf_list_to_matrix(tpt, thresh.LOD.ph =x1[i], verbose = F)
    y1[i] <- sum(!is.na(m$lod.mat))/length(m$lod.mat)
  }
  res<- rbind(res, 
              data.frame(ploidy = ploidy, n.ind = dat.sim$n.ind, LOD_thresh = x1, sparsity = y1))
}
#### BI ####
dat <- read_geno_csv(file.in = "~/repos/BI_sweetpotato.csv",
                     ploidy.p1 = 6,
                     ploidy.p2 = 6,
                     name.p1 = "Beauregard",
                     name.p2 = "Regal")
dat <- filter_missing(dat, type = "marker", filter.thres = 0.10, inter = F)
dat <- filter_missing(dat, type = "individual", filter.thres = 0.10, inter = F)
dat <- filter_individuals(dat, ind.to.remove = c("RB3", "RB6", "RB11", "RB15", "RB22", "RB30", "RB41", "RB44", "RB46", "RB50", "RB52", "RB54", "RB56", "RB58", "BR2", "BR3", "BR5", "BR7", 
                                                 "BR8", "BR9", "BR10", "BR13", "BR20", "BR21", "BR24", "BR26", "BR27", "BR28", "BR30", "BR31", "BR40", "BR41"))
s <- filter_segregation(dat, inter = FALSE)
s <- make_sequence(s, arg = "ch1")
tpt <- est_pairwise_rf(s, verbose = FALSE) ## Will remove
y1<-numeric(40)
x1<- seq(0, 40, length.out = 40)
for(i in 1:40){
  m <- rf_list_to_matrix(tpt, thresh.LOD.ph =x1[i], verbose = F)
  y1[i] <- sum(!is.na(m$lod.mat))/length(m$lod.mat)
}
res<- rbind(res, 
            data.frame(ploidy = ploidy, n.ind = dat$n.ind, LOD_thresh = x1, sparsity = y1))

#### BT ####
dat <- read_geno_csv(file.in = "~/repos/official_repos/mappoly2/misc/BT_trifida.csv",
                     ploidy.p1 = 6, ploidy.p2 = 6,
                     name.p1 = "Beauregard",
                     name.p2 = "Tanzania")
dat <- filter_missing(dat,
                      type = "marker",
                      filter.thres = 0.25, inter = F)
dat <- filter_missing(dat,
                      type = "individual",
                      filter.thres = 0.19, inter = F)
dat <- filter_individuals(dat, ind.to.remove = c("BT05.323", "BT05.324", "BT13.049", "BT13.061", "BT13.063", "BT13.135", "BT13.153", "BT13.196"))
s <- filter_segregation(dat, inter = FALSE)
s <- make_sequence(s, arg = "ch1")
s <- make_sequence(dat, s$mrk.names[sort(sample(1:1945, 500))])
tpt <- est_pairwise_rf(s, verbose = FALSE) ## Will remove
y1<-numeric(40)
x1<- seq(0, 40, length.out = 40)
for(i in 1:40){
  m <- rf_list_to_matrix(tpt, thresh.LOD.ph =x1[i], verbose = F)
  y1[i] <- sum(!is.na(m$lod.mat))/length(m$lod.mat)
}
res<- rbind(res, 
            data.frame(ploidy = ploidy, n.ind = dat$n.ind, LOD_thresh = x1, sparsity = y1))

res$n.ind <- as.factor(res$n.ind)

#### plot ####

p<-ggplot(res, aes(x=LOD_thresh, y=sparsity, group=n.ind)) +
  geom_line(aes(color=n.ind), linewidth = 1)+
  ggplot2::facet_wrap(.~ploidy) +
  scale_color_manual(values=c(viridis::viridis(6)[1], "red", viridis::viridis(6)[2:5], "orange", viridis::viridis(6)[6]))+
  xlim(0, 5)
p



#