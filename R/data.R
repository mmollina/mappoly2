#' Autotetraploid potato dataset.
#'
#' This dataset comprises 156 offspring from a cross between two tetraploid
#' potato varieties, Atlantic and B1829-5. It is genotyped with the SolCAP
#' Infinium 8303 potato array, and the genomic order of SNPs from the Solanum
#' tuberosum genome version 4.03 is included. Genotype calling was
#' performed using the ClusterCall R package. The original dataset can be found
#' in [Pereira et al., (2021)](https://doi.org/10.1038/s41437-021-00416-x)
#'
#' @format This list is of class \code{mappoly2.data} and contains the following components:
#'
#' \describe{
#' \item{ploidy.p1}{The ploidy level of parent P1 (Atlantic) is 4.}
#' \item{ploidy.p2}{The ploidy level of parent P2 (B1829-5) is 4.}
#' \item{n.ind}{The total number of individuals is 156.}
#' \item{n.mrk}{The total number of unique markers is 6511, filtered for redundancy.}
#' \item{ind.names}{A character vector of the individuals' names.}
#' \item{mrk.names}{A character vector of the markers' names.}
#' \item{name.p1}{The name of parent P1 is Atlantic.}
#' \item{name.p2}{The name of parent P2 is B1829-5.}
#' \item{dosage.p1}{A named integer vector containing the dosage in parent P1 for all \code{n.mrk} markers.}
#' \item{dosage.p2}{A named integer vector containing the dosage in parent P2 for all \code{n.mrk} markers.}
#' \item{chrom}{A named character vector indicating the chromosome each marker belongs to.}
#' \item{genome.pos}{A named numeric vector containing the physical position of the markers in the sequence.}
#' \item{ref}{A character vector of the reference allele.}
#' \item{alt}{A character vector of the alternate allele.}
#' \item{geno.dose}{A numeric matrix containing the dosage for each marker (rows) for each individual (columns).}
#' \item{redundant}{A data frame containing two character vectors, \code{kept} and \code{removed}, of markers filtered for redundancy.}
#' \item{QAQC.values}{A list containing quality assurance and quality control values with the following components:
#'   - A data frame with statistics for each marker, including `miss` (missing data rate), `chisq.pval` (chi-squared test p-value), and `read.depth` (read depth).
#'   - A data frame with statistics for each individual, including `miss` (missing data rate) and `full.sib` (indicator of non-belonging to the analyzed bi-parental cross, generated by \code{\link[mappoly]{filter_individuals}}).
#'   }
#' }
"B2721"


#' Autotetraploid alfalfa F1 dataset.
#'
#' This dataset comprises 184 offspring from a cross between two tetraploid
#' alfalfa parents, I195 and J432, which are resistant and susceptible to
#' Aphanomyces euteiches, respectively. The biparental population was genotyped
#' with the alfalfa DArTag panel described
#' in [Zhao et al., (2023)](https://www.doi.org/10.46265/genresj.EMOR6509)
#'
#' @format This list is of class \code{mappoly2.data} and contains the following components:
#'
#' \describe{
#' \item{ploidy.p1}{The ploidy level of parent P1 (I195) is 4.}
#' \item{ploidy.p2}{The ploidy level of parent P2 (J432) is 4.}
#' \item{n.ind}{The total number of individuals is 184.}
#' \item{n.mrk}{The total number of unique markers is 2795, filtered for redundancy.}
#' \item{ind.names}{A character vector of the individuals' names.}
#' \item{mrk.names}{A character vector of the markers' names.}
#' \item{name.p1}{The name of parent P1 is I195.}
#' \item{name.p2}{The name of parent P2 is J432.}
#' \item{dosage.p1}{A named integer vector containing the dosage in parent P1 for all \code{n.mrk} markers.}
#' \item{dosage.p2}{A named integer vector containing the dosage in parent P2 for all \code{n.mrk} markers.}
#' \item{chrom}{A named character vector indicating the chromosome each marker belongs to.}
#' \item{genome.pos}{A named numeric vector containing the physical position of the markers in the sequence.}
#' \item{ref}{A character vector of the reference allele.}
#' \item{alt}{A character vector of the alternate allele.}
#' \item{geno.dose}{A numeric matrix containing the dosage for each marker (rows) for each individual (columns).}
#' \item{redundant}{A data frame containing two character vectors, \code{kept} and \code{removed}, of markers filtered for redundancy.}
#' \item{QAQC.values}{A list containing quality assurance and quality control values with the following components:
#'   - A data frame with statistics for each marker, including `miss` (missing data rate), `chisq.pval` (chi-squared test p-value), and `read.depth` (read depth).
#'   - A data frame with statistics for each individual, including `miss` (missing data rate) and `full.sib` (indicator of non-belonging to the analyzed bi-parental cross, generated by \code{\link[mappoly]{filter_individuals}}).
#'   }
#' }
"alfa_f1"

#' Autotetraploid alfalfa BC dataset.
#'
#' This dataset comprises 93 offspring from a cross between two tetraploid
#' alfalfa parents, I195 and F1.85.209, the latter being derived from the cross between I195 and J432.
#' The biparental population was genotyped with the alfalfa DArTag panel described
#' in [Zhao et al., (2023)](https://www.doi.org/10.46265/genresj.EMOR6509)
#'
#' @format This list is of class \code{mappoly2.data} and contains the following components:
#'
#' \describe{
#' \item{ploidy.p1}{The ploidy level of parent P1 (I195) is 4.}
#' \item{ploidy.p2}{The ploidy level of parent P2 (F1.85.209) is 4.}
#' \item{n.ind}{The total number of individuals is 93.}
#' \item{n.mrk}{The total number of unique markers is 1923, filtered for redundancy.}
#' \item{ind.names}{A character vector of the individuals' names.}
#' \item{mrk.names}{A character vector of the markers' names.}
#' \item{name.p1}{The name of parent P1 is I195.}
#' \item{name.p2}{The name of parent P2 is F1.85.209.}
#' \item{dosage.p1}{A named integer vector containing the dosage in parent P1 for all \code{n.mrk} markers.}
#' \item{dosage.p2}{A named integer vector containing the dosage in parent P2 for all \code{n.mrk} markers.}
#' \item{chrom}{A named character vector indicating the chromosome each marker belongs to.}
#' \item{genome.pos}{A named numeric vector containing the physical position of the markers in the sequence.}
#' \item{ref}{A character vector of the reference allele.}
#' \item{alt}{A character vector of the alternate allele.}
#' \item{geno.dose}{A numeric matrix containing the dosage for each marker (rows) for each individual (columns).}
#' \item{redundant}{A data frame containing two character vectors, \code{kept} and \code{removed}, of markers filtered for redundancy.}
#' \item{QAQC.values}{A list containing quality assurance and quality control values with the following components:
#'   - A data frame with statistics for each marker, including `miss` (missing data rate), `chisq.pval` (chi-squared test p-value), and `read.depth` (read depth).
#'   - A data frame with statistics for each individual, including `miss` (missing data rate) and `full.sib` (indicator of non-belonging to the analyzed bi-parental cross, generated by \code{\link[mappoly]{filter_individuals}}).
#'   }
#' }
"alfa_bc"


