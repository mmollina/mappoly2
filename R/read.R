#'Data Input in CSV format
#'
#'Reads a comma-separated values (CSV) data file containing genetic marker data.
#'This function returns an object of class \code{mappoly2}.
#'
#' The CSV file should have rows representing markers, with the first row used
#' as a header. The first seven columns should contain the marker names, the
#' dosage in parents 1 and 2, the chromosome information (i.e. chromosome,
#' scaffold, contig, etc), the position of the marker within the sequence, the
#' alternate and reference alleles, if available. If not available, the values
#' should be filled by NA.
#' The remaining columns should contain the dosage of the full-sib population.
#' For a tetraploid example of such a file, see the \code{Examples} section.
#'
#' @param file.in a character string with the name of (or full path to)
#' the input file containing the data to be read
#'
#' @param ploidy.p1 ploidy level of parent 1
#'
#' @param ploidy.p2 ploidy level of parent 1
#'
#' @param name.p1 name of parent 1
#'
#' @param name.p2 name of parent 2
#'
#' @param filter.non.conforming if \code{TRUE} (default), data points with
#' unexpected genotypes (i.e. double reduction) are converted to 'NA'.
#' See the \code{\link[mappoly]{segreg_poly}} function for information on
#' expected classes and their respective frequencies.
#'
#' @param filter.redundant logical. If \code{TRUE} (default), removes redundant
#' markers during map construction, keeping them annotated to export to the
#' final map.
#'
#' @param verbose if \code{TRUE} (default), shows the current progress; if
#' \code{FALSE}, no output is produced
#'
#' @param x an object of the class \code{mappoly.sequence}
#'
#' @param thresh.line position of a threshold line for p values of the segregation test (default = \code{0.05/n.mrk})
#'
#' @param ... currently ignored
#'
#' @return An object of class \code{mappoly2} which contains a list with
#' the following components:
#'
#' \item{ploidy.p1}{ploidy level of the first parent}
#' \item{ploidy.p2}{ploidy level of the second parent}
#' \item{n.ind}{number individuals}
#' \item{n.mrk}{total number of markers}
#' \item{ind.names}{names or identifiers of the individuals}
#' \item{mrk.names}{names or identifiers of the genetic markers}
#' \item{name.p1}{names or identifiers of the first parent}
#' \item{name.p2}{names or identifiers of the second parent}
#' \item{dosage.p1}{the dosage for the first parent}
#' \item{dosage.p2}{the dosage for the second parent}
#' \item{chrom}{chromosome numbers for all markers}
#' \item{genome.pos}{physical positions on the genome for the genetic markers}
#' \item{seq.ref}{the reference DNA sequence data for the genetic markers}
#' \item{seq.alt}{the alternate DNA sequence data for the genetic markers}
#' \item{all.mrk.depth}{represents the depth of coverage for all genetic markers.
#' NULL when using csv input files}
#' \item{geno.dose}{a matrix containing the dosage for each marker (rows) for
#' each individual (columns).}
#' \item{kept}{a list of all non-redundant markers if \code{filter.redundant = TRUE}}
#' \item{filter.correspondence}{a list of all non-redundant markers and its equivalence
#' to the redundant ones if \code{filter.redundant = TRUE}}
#'
#' @examples
#' \donttest{
#' tempfl <- list.files(system.file('extdata', package = 'mappoly2'),
#'                      full.names = TRUE)
#' SolCAP.dose <- read_geno_csv(file.in  = tempfl,
#'                                         ploidy.p1 = 4,
#'                                         name.p1 = "Atlantic",
#'                                         name.p2 = "B1829-5")
#' print(SolCAP.dose, detailed = TRUE)
#' plot(SolCAP.dose)
#'}
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @importFrom grDevices rgb
#' @importFrom stats na.omit
#' @importFrom utils read.csv
#' @export read_geno_csv
read_geno_csv <- function(file.in,
                          ploidy.p1,
                          ploidy.p2 = ploidy.p1,
                          name.p1 = NULL,
                          name.p2 = NULL,
                          filter.non.conforming = TRUE,
                          filter.redundant = TRUE,
                          verbose = TRUE) {
  dat <- read.csv(file = file.in,
                  header = TRUE,
                  stringsAsFactors = FALSE)
  return(table_to_mappoly(dat,
                          ploidy.p1,
                          ploidy.p2,
                          name.p1,
                          name.p2,
                          filter.non.conforming,
                          filter.redundant,
                          verbose))
}

#' Conversion of data.frame to mappoly.data
#'
#' @param void internal function to be documented
#' @keywords internal
#' @export
table_to_mappoly <- function(dat,
                             ploidy.p1,
                             ploidy.p2,
                             name.p1 = NULL,
                             name.p2 = NULL,
                             filter.non.conforming = TRUE,
                             filter.redundant = TRUE,
                             verbose = TRUE) {

  # Removing markers with missing data points for parents
  dat <- dat[apply(dat[, 2:3], 1, function(x) !any(is.na(x))), ]

  # Get number of individuals
  n.ind <- ncol(dat) - 7

  # Get number of markers
  n.mrk <- nrow(dat)

  # Get marker names
  mrk.names.raw <- mrk.names <- as.character(dat[, 1, drop = TRUE])

  # Get individual's names
  ind.names <- colnames(dat)[-c(1:7)]

  # Get parent's names
  if (is.null(name.p1)) {
    name.p1 <- colnames(dat)[2]
  }
  if (is.null(name.p2)) {
    name.p2 <- colnames(dat)[3]
  }

  # Get dosage in parent P1
  dosage.p1 <- as.integer(dat[, 2, drop = TRUE])

  # Get dosage in parent P2
  dosage.p2 <- as.integer(dat[, 3, drop = TRUE])

  # Polymorphic markers
  d.p1 <- abs(abs(dosage.p1 - (ploidy.p1 / 2)) - (ploidy.p1 / 2))
  d.p2 <- abs(abs(dosage.p2 - (ploidy.p2 / 2)) - (ploidy.p2 / 2))
  id <- d.p1 + d.p2 != 0

  # Markers with parental dosages less or equal than ploidy levels
  id <- d.p1 <= ploidy.p1 & d.p2 <= ploidy.p2 & id

  # Get chromosome info
  chrom <- as.character(dat[, 4, drop = TRUE])
  chrom[is.na(chrom)] <- "NoChr"

  # Get sequence position info
  genome.pos <- as.numeric(dat[, 5, drop = TRUE])

  ref <- as.character(dat[, 6, drop = TRUE])
  alt <- as.character(dat[, 7, drop = TRUE])
  names(ref) <- names(alt) <- names(genome.pos) <- names(chrom) <- names(dosage.p2) <- names(dosage.p1) <- mrk.names

  if (verbose) {
    cat("Reading the following data:")
    txt <- list(
      paste0("    Ploidy level of ", name.p1, ":"),
      paste0("    Ploidy level of ", name.p2, ":"),
      paste0("    No. individuals:"),
      paste0("    No. markers:"),
      paste0("    No. informative markers:")
    )
    n <- sapply(txt, nchar)
    for (i in 1:length(txt)) {
      txt[[i]] <- paste(txt[[i]], paste0(rep(" ", max(n) - n[i]), collapse = ""))
    }
    cat("\n", txt[[1]], ploidy.p1)
    cat("\n", txt[[2]], ploidy.p2)
    cat("\n", txt[[3]], n.ind)
    cat("\n", txt[[4]], n.mrk)
    cat("\n ", txt[[5]], " ", sum(id), " (", round(100 * sum(id) / n.mrk, 1), "%)", sep = "")
    cat("\n     ---")
    if (any(!is.na(chrom))) {
      cat("\n     This dataset contains chromosome information.")
    }
    if (any(!is.na(genome.pos))) {
      cat("\n     This dataset contains position information.")
    }
    cat("\n     ")
  }

  # Get genotypic info
  geno.dose <- as.matrix(dat[, -c(1:7), drop = FALSE])
  dimnames(geno.dose) <- list(mrk.names, ind.names)

  # Replacing non-conforming values with NA
  geno.dose[geno.dose > ploidy.p1/2 + ploidy.p2/2] <- NA

  geno.dose <- geno.dose[id, , drop = FALSE]
  res <- structure(
    list(data = list(
      ploidy.p1 = ploidy.p1,
      ploidy.p2 = ploidy.p2,
      n.ind = n.ind,
      n.mrk = sum(id),
      ind.names = ind.names,
      mrk.names = mrk.names[id],
      name.p1 = name.p1,
      name.p2 = name.p2,
      dosage.p1 = dosage.p1[id],
      dosage.p2 = dosage.p2[id],
      chrom = chrom[id],
      genome.pos = genome.pos[id],
      ref = ref,
      alt = alt,
      all.mrk.depth = NULL,
      geno.dose = geno.dose,
      redundant = NULL)),
    class = c("mappoly2", "data")
  )

  # Screening non-conforming markers
  if (filter.non.conforming) {
    if (verbose) cat(" -->  Filtering non-conforming markers.\n     ")
    res$data <- mappoly2:::filter_non_conforming_classes(res$data)
  }
  # Screening redundant markers
  if(filter.redundant){
    if (verbose) cat(" -->  Filtering markers with redundant information.\n     ")
    redundant <- mappoly2:::filter_redundant(res)
    if(all(is.na(redundant))) res$data$redundant <- NA
    else{
      res$data <- mappoly2:::subset_data(res$data, select.mrk = setdiff(res$data$mrk.names, redundant$removed))
    }
    res$data$redundant <- redundant
  }
  res$data$QAQC.values <- .setQAQC(id.mrk = res$data$mrk.names,
                                   id.ind = res$data$ind.names,
                                   miss.mrk = apply(res$data$geno.dose, 1, function(x) sum(is.na(x)))/res$data$n.ind,
                                   miss.ind = apply(res$data$geno.dose, 2, function(x) sum(is.na(x)))/res$data$n.mrk,
                                   chisq.pval = suppressWarnings(mappoly2:::mappoly_chisq_test(res$data)))
  if(verbose) cat("----------------------------------\n")
  return(res)
}
