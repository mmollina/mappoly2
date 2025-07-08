#' Autotetraploid Potato Dataset
#'
#' This dataset contains 156 offspring derived from a cross between two tetraploid
#' potato varieties, Atlantic and B1829-5. It was genotyped using the SolCAP Infinium
#' 8303 potato array, and the genomic order of SNPs follows the Solanum tuberosum
#' genome version 4.03. Genotype calling was performed with the ClusterCall R package.
#' The original dataset is available in [Pereira et al. (2021)](https://doi.org/10.1038/s41437-021-00416-x).
#'
#' @format A list of class \code{mappoly2.data} containing the following components:
#' \describe{
#'   \item{ploidy.p1}{The ploidy level of parent P1 (Atlantic): 4.}
#'   \item{ploidy.p2}{The ploidy level of parent P2 (B1829-5): 4.}
#'   \item{n.ind}{The total number of individuals: 156.}
#'   \item{n.mrk}{The number of unique markers after filtering for redundancy: 6511.}
#'   \item{ind.names}{Character vector containing the names of the individuals.}
#'   \item{mrk.names}{Character vector containing the names of the markers.}
#'   \item{name.p1}{The name of parent P1: Atlantic.}
#'   \item{name.p2}{The name of parent P2: B1829-5.}
#'   \item{dosage.p1}{Named integer vector with the dosage for parent P1 across \code{n.mrk} markers.}
#'   \item{dosage.p2}{Named integer vector with the dosage for parent P2 across \code{n.mrk} markers.}
#'   \item{chrom}{Named character vector indicating the chromosome each marker belongs to.}
#'   \item{genome.pos}{Named numeric vector with the physical positions of the markers.}
#'   \item{ref}{Character vector of reference alleles.}
#'   \item{alt}{Character vector of alternate alleles.}
#'   \item{geno.dose}{Numeric matrix containing dosage information for each marker (rows) and individual (columns).}
#'   \item{redundant}{Data frame listing \code{kept} and \code{removed} markers filtered for redundancy.}
#'   \item{QAQC.values}{A list containing quality assurance and quality control statistics:
#'     \describe{
#'       \item{Marker Statistics}{Data frame with statistics such as:
#'         \code{miss} (missing data rate),
#'         \code{chisq.pval} (chi-squared test p-value),
#'         and \code{read.depth} (read depth).}
#'       \item{Individual Statistics}{Data frame with statistics for each individual:
#'         \code{miss} (missing data rate),
#'         and \code{full.sib} (indicator of misclassified individuals not belonging
#'         to the cross, generated using \code{\link[mappoly]{filter_individuals}}).}
#'     }
#'   }
#' }
"potato"

#' Autotetraploid Alfalfa F1 Dataset
#'
#' This dataset contains 184 offspring from a cross between two tetraploid
#' alfalfa parents, I195 (resistant to Aphanomyces euteiches) and J432 (susceptible).
#' The population was genotyped using the alfalfa DArTag panel described in
#' [Zhao et al. (2023)](https://www.doi.org/10.46265/genresj.EMOR6509).
#'
#' @format A list of class \code{mappoly2.data} with the following components:
#' (Structure is identical to the potato dataset above, with specifics adapted.)
"alfalfa_f1"

#' Autotetraploid Alfalfa BC Dataset
#'
#' This dataset includes 93 offspring from a backcross between tetraploid alfalfa
#' parents I195 and F1.85.209 (a hybrid derived from I195 and J432). Genotyping was
#' performed using the alfalfa DArTag panel as described in
#' [Zhao et al. (2023)](https://www.doi.org/10.46265/genresj.EMOR6509).
#'
#' @format A list of class \code{mappoly2.data} with the following components:
#' (Structure is identical to the potato dataset above, with specifics adapted.)
"alfalfa_bc"

#' Diploid Apple F1 Dataset – 4-Chromosome Subset
#'
#' This dataset contains SNP dosage calls for **318 F₁ seedlings** obtained from the
#' diploid apple cross *‘Jonathan’ × ‘Golden Delicious’*.
#' The 610 bi-allelic SNPs included here belong to **four chromosomes
#' (Chr 2, 3, 8 and 9)** and were generated with restriction-site–associated DNA
#' sequencing (RAD-seq) as described in Sun *et al.* (2015).
#' The subset is provided for tutorials and vignette examples in **\pkg{mappoly2}**.
#'
#' ## Details
#' * **Parents:** ‘Jonathan’ (maternal) and ‘Golden Delicious’ (paternal)
#' * **Population size:** 318 individuals
#' * **Markers:** 610 SNPs
#'   * Chr 2 – 134; Chr 3 – 160; Chr 8 – 149; Chr 9 – 167
#' * **Ploidy:** Diploid (2 × = 34)
#'
#' ## Source
#' Sun R, Chang Y, Yang F, *et al.* (2015)
#' *A dense SNP genetic map constructed using restriction site-associated DNA
#' sequencing enables detection of QTLs controlling apple fruit quality.*
#' **BMC Genomics 16**:747. <doi:10.1186/s12864-015-1946-x>
#'
#' @format A list of class \code{mappoly2.data} with the following
#' components.
#' *(Structure is identical to the \code{alfalfa_f1} dataset, with these
#' apple-specific values.)*
#'
#' @keywords datasets
"apple"

