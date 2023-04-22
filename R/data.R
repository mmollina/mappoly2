#' Autotetraploid potato dataset.
#'
#' This dataset comprises 156 offsprings of a cross between two tetraploid
#' potato varieties, Atlantic and B1829-5. It is genotyped with the SolCAP
#' Infinium 8303 potato array, and the genomic order of SNPs from the Solanum
#' tuberosum genome version 4.03 is included. The genotype calling was
#' performed using the ClusterCall R package. The original dataset can be found
#' in [Pereira et al., (2021)](https://doi.org/10.1038/s41437-021-00416-x)
#'
#' @format This list is of class \code{mappoly2.data} and contains the following components:
#'
#' \describe{
#' \item{ploidy.p1}{The ploidy level of parent P1 is 4.}
#' \item{ploidy.p2}{The ploidy level of parent P2 is 4.}
#' \item{n.ind}{The total number of individuals is 156.}
#' \item{n.mrk}{The total number of markers is 6511, filtered for redundancy.}
#' \item{ind.names}{A character vector of the names of the individuals.}
#' \item{mrk.names}{A character vector of the names of the markers.}
#' \item{name.p1}{The name of parent P1 is Atlantic.}
#' \item{name.p2}{The name of parent P2 is B1829-5.}
#' \item{dosage.p1}{A named integer vector containing the dosage in parent P1 for all \code{n.mrk} markers.}
#' \item{dosage.p2}{A named integer vector containing the dosage in parent P2 for all \code{n.mrk} markers.}
#' \item{chrom}{A named character vector indicating the chromosome each marker belongs to.}
#' \item{genome.pos}{A named numeric vector containing the physical position of the markers in the sequence.}
#' \item{ref}{A character vector of the reference allele.}
#' \item{alt}{A character vector of the alternate allele.}
#' \item{prob.thres}{A null value.}
#' \item{geno.dose}{A numeric matrix containing the dosage for each marker (rows) for each individual (columns).}
#' \item{chisq.pval}{A named numeric vector containing p-values for all markers associated with the chi-square test for the expected segregation patterns under Mendelian segregation.}
#' \item{redundant}{A data frame containing two character vectors, \code{kept} and \code{removed}, of markers filtered for redundancy.}
#' }
"B2721"
