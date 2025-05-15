# MAPpoly2 - Under development package

MAPpoly2 is an R package designed to build maps in interconnected full-sib 
autopolyploid families. It has been enhanced for performence, user-friendliness 
and accessibility. This version, developed for potential integration with R Shiny, 
aims to provide a user-intuitive interface for genetic mapping in polyploids. 
It can handle ploidy levels of 2, 4, and 6, including any combination of these.

One of the key improvements in MAPpoly2 is its enhanced performance, largely due 
to the implementation of computationally intensive codes primarily in C++. This 
enables efficient handling of large datasets. Additionally, the package implements 
the construction of individual maps for each parent using a Hidden Markov Model (HMM), 
significantly speeding up the map construction process. These individual maps can 
then be merged, and a joint map is recomputed to include any remaining markers.

## Installation:

## From GitHub 

You can install the development version from Git Hub. Within R, you need to 
install `devtools`:

```R
install.packages("devtools")
```

If you are using Windows, please install the latest recommended 
version of [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

To install MAPpoly from 2025_updates Git Hub branch, use

```R
devtools::install_github("mmollina/mappoly2@2025_updates", dependencies=TRUE)
```

## Tutorial

- [Building an integrated genetic linkage map of autotetraploid alfalfa populations using the MAPpoly2](https://rpubs.com/mmollin/tutorial_mappoly2)

- [Code Only - Alfalfa populations using the MAPpoly2](https://github.com/mmollina/mappoly2_vignettes/blob/main/mappoly2_alfalfa.R) 

- [Constructing Multi-Family Genetic Maps with MAPpoly2: A Simulation Example](https://rpubs.com/mmollin/multi_family_simulation)

- [Supplementary Slides](https://github.com/mmollina/mappoly2_vignettes/blob/main/Updates-Introducing_MAPpoly2-and_updates_QTLpoly-2024-workshop.pdf)


# Task List

## Testing & Debugging
- Develop a comprehensive test script.
- Identify and resolve existing bugs.

## Algorithm Development
- ~~Implement an algorithm similar to `mappoly::est_rf_hmm_sequential`~~.
- Design an algorithm for reconstructing offspring haplotypes, 
detailing crossover points and identifying homologs involved 
in the exchange.

## Data Filtering
- Split the existing data filtering function into specific filters for:
  - ~~Individual data~~
  - ~~Marker data~~
  - ~~Segregation data~~
  - ~~Read depth~~
- Filter individuals:
  - Remove “x” in parents on the plot.
  - Check missing rate and apply if necessary.
- Use a marker and recompute χ² when needed (if not filtered by χ²).
- List markers removed due to filters.
- Include redundant markers back.

## Phase Estimation & Adjustment
- Implement a joint phase estimation method for the entire population, 
incorporating a two-point linkage and algorithm.
- Make phase information available throughout the process.
- Develop a function to split and rephase data as needed.

## Data Integration
- Enable support for reading **DarTAG** data.
- ~~Add functionality for reading multiparental population data~~.
- Implement plot method for multiparental dataset
- Create a function for reading diploid data.

## Population Types
- Extend support to include inbred-based diploid populations.

## Map Augmentation
- Develop functionality for augmenting maps for single-parent datasets.
- Support the integration of phase data from one order into another.
- Ensure "pi/pi" markers are included when augmenting maps.

## Visualization
- Implement multiple plotting and printing methods to enhance data visualization.
- ~~Plot MDS (Multidimensional Scaling) to visualize marker relationships~~.

## Documentation
- Create detailed documentation outlining the structure and components of the package objects.

## External Package Connectivity
- Establish connections with **AlphaSimR** for simulation integration.
- Enable compatibility with **R/qtl** for QTL analysis.

## Data Validation
- Verify segregation p-values for accuracy.

## Marker Grouping
- Provide an easy way to group markers.

# Acknowledgment

This package has been developed as part of the project [AFRI-Grant: A Genetics-Based Data Analysis System for Breeders in Polyploid Breeding Programs](https://portal.nifa.usda.gov/web/crisprojectpages/1027948-a-genetics-based-data-analysis-system-for-breeders-in-polyploid-breeding-programs.html) and  [SCRI-Grant: Tools for polyploids](https://www.polyploids.org/), funded by USDA NIFA.


<div class="horizontalgap" style="width:5px">
     <a id="USDA-NIFA" href="https://portal.nifa.usda.gov/web/crisprojectpages/1027948-a-genetics-based-data-analysis-system-for-breeders-in-polyploid-breeding-programs.html"><img src="nifa-color-lockup.png" width="650" alt=""/></a> 
      <a id="NCSU" href="https://www.ncsu.edu/"><img src="https://brand.ncsu.edu/assets/logos/ncstate-brick-2x2-red.png" width="150" alt=""/></a>
    <span class="stretch"></span>
</div>


---
<sub>NC State University promotes equal opportunity and prohibits discrimination and harassment based upon one’s age, color, disability, gender identity, genetic information, national origin, race, religion, sex (including pregnancy), sexual orientation and veteran status.</sub>



