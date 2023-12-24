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

If you are using Windows, please install the the latest recommended version of [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

To install MAPpoly from Git Hub use

```R
devtools::install_github("mmollina/mappoly", dependencies=TRUE)
```

# Acknowledgment

This package has been developed as part of the [Genomic Tools for Sweetpotato Improvement project](https://sweetpotatogenomics.cals.ncsu.edu/) (GT4SP) and [SweetGAINS](https://cgspace.cgiar.org/handle/10568/106838), both funded by [Bill & Melinda Gates Foundation](https://www.gatesfoundation.org/). Its continuous improvement is made possible by the project [AFRI-Grant: A Genetics-Based Data Analysis System for Breeders in Polyploid Breeding Programs](https://portal.nifa.usda.gov/web/crisprojectpages/1027948-a-genetics-based-data-analysis-system-for-breeders-in-polyploid-breeding-programs.html) and  [SCRI-Grant: Tools for polyploids](https://www.polyploids.org/), funded by USDA NIFA.

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



