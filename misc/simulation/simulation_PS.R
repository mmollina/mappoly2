sim_map  <- function(ploidy = 4,                       # ploidy level
                     map.length = 100,                 # distance in cM
                     n.ind = 200 ,                     # number of individuals
                     n.mrk = 200,                      # number of markers
                     n.chr = 3,                        # number of chromosomes
                     perc.random.missing = 10,        # percentage on random missing data
                     prob.dose = c(0.225,0.225,0.1,0.225,0.225),
                     # dose frequency for each parent [Pr(dose=0), ..., Pr(dose=ploidy)]
                     seed.for.config = 40,             # seed for simulating the linkage phase configuration
                     seed.for.pop = 1,                 # seed for simulating recombinations
                     prefPairing  =  0,                # preferential pairing
                     quadrivalents  =  0,              # quadrivalent formation rate
                     centromere  =  20,                # long arm 4 times as long as the short arm
                     natural.pairing  =  0)            # probabilities of bivalents and quadrivalents are determined by the “quadrivalents” argument
{
  genome.pos <- chrom <- w <- NULL
  for(i in 1:n.chr)
  {
    ph.temp<-mappoly::sim_homologous(ploidy = ploidy,
                                     n.mrk=n.mrk,
                                     max.d=ploidy,    #max dosage number
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
  P<-as.matrix(G[,c(1:2)])
  colnames(P) <- c("P1", "P2")
  G<-G[,-c(1:2)]
  colnames(G)<-paste0("Ind_", 1:n.ind)
  ref <- sample(c("A", "T", "C", "G"), replace = T, size = nrow(P))
  alt <- character(length(ref))
  for(i in 1:length(ref))
    alt[i] <- sample(setdiff(c("A", "T", "C", "G"), ref[i]), size = 1)
  if(perc.random.missing > 0){
    n <- ceiling(length(G) * perc.random.missing/100)
    G[sample(1:length(G),n)] <- NA
  }
  DF <- data.frame(snp_id = rownames(P),
                   P,
                   chrom = chrom,
                   genome_pos = genome.pos,
                   ref = ref,
                   alt = alt,
                   G)
  dat <- mappoly2::table_to_mappoly(dat = DF, ploidy.p1 = ploidy, ploidy.p2 = ploidy, verbose = FALSE)
  return(dat)
}
