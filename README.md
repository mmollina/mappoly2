# MAPpoly2

MAPpoly2 is an R package designed to build maps in interconnected full-sib autopolyploid families. It has been enhanced for user-friendliness and accessibility. This version, developed for potential integration with R Shiny, aims to provide a user-intuitive interface for genetic mapping in polyploids. It can handle ploidy levels of 2, 4, and 6, including any combination of these.

One of the key improvements in MAPpoly2 is its enhanced performance, largely due to the implementation of computationally intensive codes primarily in C++. This enables efficient handling of large datasets. Additionally, the package facilitates the construction of individual maps for each parent using a Hidden Markov Model (HMM), significantly speeding up the map construction process. These individual maps can then be merged, and a joint map is recomputed to include any remaining markers.

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
[Constructing Multi-Family Genetic Maps with MAPpoly2: A Simulation Example](https://rpubs.com/mmollin/multi_family_simulation)

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



