rm(list = ls())
require(mappoly2)
P1 <- "AC46"
P2 <- "AC78"
x <- read_geno_csv(file.in = paste0("~/repos/collaborations/acacia/full_sib_data/mappoly2_", P1, "-", P2,".csv"),
                     ploidy.p1 = 2,
                     ploidy.p2 = 2,
                     name.p1 = P1,
                     name.p2 = P2)
print(x, detailed = TRUE)
plot(x, what = "raw")
x <- filter_data(x, mrk.thresh = .08, ind.thresh = .7)
plot(x, what = "screened")
x <- set_initial_sequence(x, c("ch1", "ch13"))
plot(x, what = "initiated")
x <- pairwise_rf(x, ncpus = 7)
plot(x, what = "pairwise", fact = 1)
x <- init_rf_filter(x, 5, 5, .1, c(0.02, 0.98))
plot(x, what = "pairwise", fact = 1, ord = x$initial.screened.rf$mrk.names)
x <- group(x, expected.groups = 2, inter = FALSE)
plot(x, what = "group")
x <- set_working_sequence(x, ch = list("ch1", "ch13"))
x <- genome_order(x)
plot(x, what = "pairwise_genome")
x <- mds(x)
plot(x, what = "pairwise_mds")
x <- pairwise_phasing_per_group(x,
                                gr = 1,
                                type = "genome",
                                thresh.LOD.ph = 10,
                                thresh.LOD.rf = 10)
x <- mapping_per_group(x, type = "genome", gr = 1, error = 0.05, verbose = TRUE)
x <- pairwise_phasing_per_group(x,
                                gr = 2,
                                type = "genome",
                                thresh.LOD.ph = 10,
                                thresh.LOD.rf = 10)
x <- mapping_per_group(x, type = "genome", gr = 2, error = 0.05, verbose = TRUE)
x <- pairwise_phasing_per_group(x,
                                gr = 1,
                                type = "mds",
                                thresh.LOD.ph = 10,
                                thresh.LOD.rf = 10)
x <- mapping_per_group(x, type = "mds", gr = 1, error = 0.05, verbose = TRUE)
x <- pairwise_phasing_per_group(x,
                                gr = 2,
                                type = "mds",
                                thresh.LOD.ph = 10,
                                thresh.LOD.rf = 10)
x <- mapping_per_group(x, type = "mds", gr = 2, error = 0.05, verbose = TRUE)
plot_map(x, gr = 1, type = "mds")
plot_map(x, gr = 1, type = "genome")
plot_map(x, gr = 2, type = "mds")
plot_map(x, gr = 2, type = "genome")
x <- calc_haplotypes_per_group(x, 1)
x <- calc_haplotypes_per_group(x, 2)

#image(t(as.matrix(x$working.sequences$Seq_1$order$genome$hmm.map$haploprob[,-c(1:2)])))
#image(t(as.matrix(x$working.sequences$Seq_1$order$mds$hmm.map$haploprob[97:98,-c(1:2)])))
#image(t(as.matrix(x$working.sequences$Seq_1$order$genome$hmm.map$haploprob[97:98,-c(1:2)])))

x <- augment_phased_map_per_group(x, gr = 1, type = "mds")
x <- augment_phased_map_per_group(x, gr = 1, type = "genome")
x <- augment_phased_map_per_group(x, gr = 2, type = "mds")
x <- augment_phased_map_per_group(x, gr = 2, type = "genome")

plot_map(x, gr = 1, type = "mds")
plot_map(x, gr = 1, type = "genome")
plot_map(x, gr = 2, type = "mds")
plot_map(x, gr = 2, type = "genome")

M1 <- as.matrix(x$working.sequences$Seq_1$order$genome$hmm.map$haploprob)
M1[1:10, 1:10]
M1[,2] <- rep(1:4, nrow(M1)/4)
M1[1:10, 1:10]
for(i in 1:nrow(M1)){
  z <- M1[i,]
  M1[i,-c(1:2)] <- ifelse(z[-c(1:2)] > 0.9, z[2], NA)
}
M1 <- M1[,-c(1,2)]
par(bg = 1)
scales::show_col(ggsci::pal_tron("legacy")(4), labels = FALSE)
text(x = c(0.5, 1.5, 0.5, 1.5), y = c(-.5,-.5,-1.5, -1.5),
     labels = c("P1_A", "P1_a", "P2_A", "P2_a"), cex = 3)
image(t(M1), col = ggsci::pal_tron()(4))
image(t(M1[1:200,]), col = ggsci::pal_tron()(4))
