# MAPpoly2

MAPpoly2 is an R package that builds upon the successful foundation of MAPpoly, 
specifically tailored for genetic mapping in polyploids. It represents an advancement,
especially in managing large datasets and in facilitating interactive applications. 

As an early-stage developmental version, MAPpoly2 has been 
designed to enhance user-friendliness and accessibility. This version has been developed
for potential integration with R Shiny, aiming to offer a dynamic and user-intuitive
interface for genetic mapping in polyploids. 

One of the key improvements in MAPpoly2 is its enhanced performance, largely attributable to the implementation of computationally intensive codes primarily in C++. This allows for efficient handling of large datasets. Additionally, the package enables the construction of individual maps for each parent using a Hidden Markov Model (HMM), which significantly speed up the map construction process. These individual maps can be subsequently merged, and a joint map is recomputed to include any remaining markers.

Moreover, MAPpoly2 introduces a more streamlined mapping process. This improvement is particularly significant as it lays the groundwork for future integration with Shiny.

## Main Functions:
- `add_marker`: Add markers to a pre-mapped sequence
- `augment_phased_map`: Augment a phased map with additional information
- `drop_marker`: Remove markers from a sequence
- `pairwise_phasing`: Perform pairwise phasing on sequences
- `filter_data`, `filter_individuals`: Functions for data filtering
- `genome_order`, `group`: Functions for ordering genome sequences and grouping
- `make_sequence`, `mapping`: Create and map genetic sequences
- `plot.mappoly2.data`, `plot_genome_vs_map`, `plot_map`: Visualization functions
- `read_geno_csv`: Read genotype data from CSV
- `rev_map`, `rf_filter`: Reverse mapping and filter recombination fractions
- Other utility functions like `subset.mappoly2.data`, `plot_map_list`


## Installation:

## From GitHub 

You can install the development version from Git Hub. Within R, you need to install `devtools`:

```R
install.packages("devtools")
```

If you are using Windows, please install the latest recommended version of [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

To install MAPpoly from Git Hub, use

```R
devtools::install_github("mmollina/mappoly", dependencies=TRUE)
```


## Tutorial

[Building an integrated genetic linkage map of autotetraploid alfalfa populations using the MAPpoly2](https://rpubs.com/mmollin/tutorial_mappoly2)

## Test code

```R
# Clear the workspace
rm(list = ls())
require(mappoly2)

# Set up ploidy and parent information
ploidy.vec <- c(4, 2, 4, 2, 4, 4) # Three parents
names(ploidy.vec) <- c("P1", "P2", "P3", "P4", "P5", "P6")

# Define parent pairs for crosses
parents.mat <- matrix(c("P1", "P2", 
                        "P1", "P3", 
                        "P2", "P2", 
                        "P3", "P4", 
                        "P4", "P5", 
                        "P5", "P6", 
                        "P5", "P2"),
                      ncol = 2, byrow = TRUE)
n.mrk <- rep(500, length(ploidy.vec)) # Number of markers
alleles <- list(P1 = 0:1, P2 = 0:1, P3 = 0:1, P4 = 0:1, P5 = 0:1, P6 = 0:1) # Alleles
n.ind <- rep(200, nrow(parents.mat)) # Number of individuals
n.chrom <- 8 # Number of chromosomes
map.length <- round(runif(n.chrom, 60, 100)) # Length of maps

# Simulate multiple crosses for each chromosome
x <- vector("list", n.chrom)
for (i in 1:n.chrom) {
  cat("chrom: ", i, "\n")
  x[[i]] <- mappolymp:::simulate_multiple_crosses(ploidy.vec, parents.mat, n.ind, n.mrk, alleles, map.length[i])
}

# Convert simulation output to mappoly2 format
source(file = "https://raw.githubusercontent.com/mmollina/MAPpoly_vignettes/master/misc_scripts/mappoly_mp_simulation_to_mappoly2.R")
system("rm ~/repos/official_repos/misc_test_mp2/multipop/*")
mappoly_mp_simulation_to_mappoly2(x, file_path = "~/repos/official_repos/misc_test_mp2/multipop")

# Rename parents for visualization
pm <- parents.mat
pm[parents.mat == "P1"] <- "Santa"
pm[parents.mat == "P2"] <- "Claus"
pm[parents.mat == "P3"] <- "Eastern"
pm[parents.mat == "P4"] <- "Bunny"
pm[parents.mat == "P5"] <- "Tooth"
pm[parents.mat == "P6"] <- "Fairy"

# Process and visualize the maps for each population
MAPs <- vector("list", length(list))
for (i in 1:nrow(parents.mat)) {
  p1 <- pm[i,1]
  p2 <- pm[i,2]
  fl <- paste0("~/repos/official_repos/misc_test_mp2/multipop/",parents.mat[i,1],"x",parents.mat[i,2],".csv")
  dat <- read_geno_csv(file.in = fl,
                       ploidy.p1 = ploidy.vec[parents.mat[i,1]],
                       ploidy.p2 = ploidy.vec[parents.mat[i,2]],
                       name.p1 = p1, name.p2 = p2)
  dat <- filter_data(dat)
  dat <- pairwise_rf(dat, mrk.scope = "per.chrom")
  dat <- rf_filter(dat, probs = c(0, 1))
  g <- group(x = dat, expected.groups = n.chrom, comp.mat = TRUE, inter = FALSE)
  s <- make_sequence(g)
  s <- order_sequence(s, type = "genome")
  s <- pairwise_phasing(s, type = "genome", parent = "p1p2")
  s <- mapping(s, type = "genome", parent = "p1p2", ncpus = n.chrom)
  s <- calc_haplotypes(s, type = "genome", ncpus = n.chrom)
  
  # Append the processed maps to the MAPs list
  MAPs[[i]] <- s
}

# Visualization of maps and haplotypes
plot_map_list(MAPs[[1]], type = "genome", col = mappoly::mp_pallet3(n.chrom))
plot_haplotypes(MAPs[[1]], type = "genome", lg = 1)
plot_map(MAPs[[1]], type = "genome")
map_summary(MAPs[[1]], type = "genome")
plot_multi_map(MAPs)

# Prepare and visualize integrated maps
x <- prepare_to_integrate(x = MAPs)
plot(x)
x <- estimate_consensus_map(x, ncpus = n.chrom)
plot(x)
plot(x, only.consensus = TRUE, col = mappoly::mp_pallet3(n.chrom))

# Calculate and visualize consensus haplotypes
system.time(x <- calc_consensus_haplo(x, ncpus = n.chrom))
plot_consensus_haplo(x, lg = 1, ind = "Ind_P1xP2_1")
plot_consensus_haplo(x, lg = 1, ind = "Ind_P2xP2_1")
```

# Acknowledgment

This package has been developed as part of the project [AFRI-Grant: A Genetics-Based Data Analysis System for Breeders in Polyploid Breeding Programs](https://portal.nifa.usda.gov/web/crisprojectpages/1027948-a-genetics-based-data-analysis-system-for-breeders-in-polyploid-breeding-programs.html) and  [SCRI-Grant: Tools for polyploids](https://www.polyploids.org/), funded by USDA NIFA.

<div class="horizontalgap" style="width:5px">
    <a id="NCSU" href="https://www.ncsu.edu/"><img src="https://brand.ncsu.edu/assets/logos/ncstate-brick-2x2-red.png" width="150" alt=""/></a>
    <a id="BMGF" href="https://www.gatesfoundation.org/"><img src="https://fsm-alliance.org/wp-content/uploads/gates-logo-bda5cc0866e8e37eccab4ac502b916c1-copy.png" width="150" alt=""/></a>
    <a id="GT4SP" href="https://sweetpotatogenomics.cals.ncsu.edu/"><img src="http://www.sweetpotatoknowledge.org/wp-content/uploads/2016/02/GT4SP-logo-e1456736272456.png" width="70" alt=""/></a>
    <a id="sweetgains" href="https://cgspace.cgiar.org/handle/10568/106838"><img src="https://cipotato.org/wp-content/uploads/2020/06/SweetGains-sin-fondo-1-350x230.png" width="150" alt=""/></a>
    <a id="PolyploidTools" href="https://www.polyploids.org/"><img src="https://www.polyploids.org/sites/default/files/inline-images/Project%20Logo-transparent.png" width="180" alt=""/></a>    
     <a id="USDA-NIFA" href="https://portal.nifa.usda.gov/web/crisprojectpages/1027948-a-genetics-based-data-analysis-system-for-breeders-in-polyploid-breeding-programs.html"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/USDA_NIFA_Twitter_Logo.jpg" width="100" alt=""/></a>  
    <span class="stretch"></span>
</div>

---
<sub>NC State University promotes equal opportunity and prohibits discrimination and harassment based upon one’s age, color, disability, gender identity, genetic information, national origin, race, religion, sex (including pregnancy), sexual orientation and veteran status.</sub>



