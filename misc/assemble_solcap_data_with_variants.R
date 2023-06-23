a <- read.csv("~/repos/official_repos/mappoly2/misc/B2721_dose.csv")
a[1:10,1:10]
b <-  read.csv("~/repos/official_repos/mappoly2/misc/solcal_ref_alt.csv")
b[1:10, ]
d<-stringr::str_split_fixed(b[,3], "\\[|\\/|\\]", n = 4)[,2:3]
rownames(d) <- b[,2]
a$sequence[is.na(a$sequence)]<-0
x<-data.frame(snp_id = a$snp_name,
              P1 = a$P1,
              P2 = a$P2,
              chrom = paste0("Chr_", sprintf("%02d", a$sequence)),
              genome4_pos = a$sequence_position,
              REF = d[a$snp_name,1],
              ALT = d[a$snp_name,2],
              a[,-c(1:5)])
write.csv(x, "~/repos/official_repos/mappoly2/inst/extdata/B2721.csv", row.names = F)
