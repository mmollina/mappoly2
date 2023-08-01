require(tidyverse)
dose.Dart <- read.csv("~/repos/collaborations/sweetpotato-breeding-insight/corrected_data/DSp22-7577_Allele_Dose_Report_updateID.csv",
                      skip = 7, row.names = 1)
counts.Dart <- read.csv("~/repos/collaborations/sweetpotato-breeding-insight/corrected_data/DSp22-7577_Allele_match_counts_collapsed_updateID.csv",
                        skip = 7, row.names = 1)
f1 <- dose.Dart[,colnames(dose.Dart)[str_detect(colnames(dose.Dart), "BR") | str_detect(colnames(dose.Dart), "RB")]]
p1<- dose.Dart[,colnames(dose.Dart)[str_detect(colnames(dose.Dart), "BEAUREGARD")]]
p2<- dose.Dart[,colnames(dose.Dart)[str_detect(colnames(dose.Dart), "REGAL")]]
#### Using markers where the two replicates of Beauregard had the same dosage calling
mrk.id <- which(apply(p1, 1, function(x) length(unique(x))==1))
#### Gathering genome positiomn
genome.pos <- as.numeric(sapply(strsplit(names(mrk.id), split = "Chr|_"), function(x) x[3]))
chrom <- as.numeric(sapply(strsplit(names(mrk.id), split = "Chr|_"), function(x) x[2]))

# Function to find different letter
find_SNP <- function(s){
  for(i in 1:nchar(s[1])){
    if(substr(s[1], i, i) != substr(s[2], i, i)){
      return(c(ref = substr(s[1], i, i), alt = substr(s[2], i, i)))
    }
  }
  return(NULL)
}
varinats <- matrix(NA, length(mrk.id), 2, dimnames = list(names(mrk.id), c("ref", "alt")))
for(i in names(mrk.id)){
  varinats[i,] <- find_SNP(s = counts.Dart[counts.Dart$CloneID %in% i, "AlleleSequence"])
}
#### Data frame form MAPpoly
DF <- cbind(snp_id = names(mrk.id),
            P1 = p1[mrk.id,1],
            P2 = p2[mrk.id],
            chrom = chrom,
            genome_pos = genome.pos,
            ref = varinats[names(mrk.id),"ref"],
            alt = varinats[names(mrk.id),"alt"],
            f1[mrk.id,])
write.csv(DF, file = "~/repos/BI_sweetpotato.csv", quote = FALSE, row.names = FALSE)
