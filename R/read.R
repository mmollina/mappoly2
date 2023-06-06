#'Data Input in CSV format
#'
#'Reads a comma-separated values (CSV) data file containing genetic marker data.
#'This function returns an object of class \code{mappoly.data}.
#'
#' The CSV file should have rows representing markers, with the first row used
#' as a header. The first seven columns should contain the marker names, the
#' dosage in parents 1 and 2, the chromosome information (i.e. chromosome,
#' scaffold, contig, etc), the position of the marker within the sequence, the
#' alternate and reference alleles, if available.
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
#' @return An object of class \code{mappoly.data} which contains a list with
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
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu} and
#' Gabriel Gesteira, \email{gdesiqu@ncsu.edu}
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
  # swap.parents <- FALSE
  # if(ploidy.p1 > ploidy.p2){
  #   swap.parents <- TRUE
  #   ## swap ploidy levels if p1 > p2
  #   temp <- ploidy.p2
  #   ploidy.p2 <- ploidy.p1
  #   ploidy.p1 <- temp
  #   ## swap names if p1 > p2
  #   temp <- name.p2
  #   name.p2 <- name.p1
  #   name.p1 <- temp
  #   ## swap dosages
  #   dat[,2:3] <- dat[,3:2]
  #   names(dat)[2:3] <- names(dat)[3:2]
  # }

  # Removing markers with missing data points for parents
  dat <- dat[apply(dat[, 2:3], 1, function(x) !any(is.na(x))), ]

  # Get number of individuals
  n.ind <- ncol(dat) - 7

  # Get number of markers
  n.mrk <- nrow(dat)

  # Get marker names
  mrk.names <- as.character(dat[, 1, drop = TRUE])

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
    list(
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
      prob.thres = NULL,
      geno.dose = geno.dose,
      chisq.pval = NULL,
      redundant = NULL#,
      #swap.parents = swap.parents
    ),
    class = "mappoly2.data"
  )
  # Computing chi-square p.values
  res <- suppressWarnings(mappoly_chisq_test(res))

  # Screening non-conforming markers
  if (filter.non.conforming) {
    if (verbose) cat(" -->  Filtering non-conforming markers.\n     ")
    res <- filter_non_conforming_classes(res)
  }
  # Screening redundant markers
  if(filter.redundant){
    if (verbose) cat(" -->  Filtering redundant markers.\n     ")
    s <- make_sequence(res, arg = "all")
    sf <- filter_redundant(s)
    res$redundant <- sf$redundant
    res <- subset_data(res, select.mrk = setdiff(res$mrk.names, sf$redundant$removed))
  }
  if(verbose) cat("----------------------------------\n")
  return(res)
}
#' @rdname read_geno_csv
#' @export
print.mappoly2.data <- function(x, ...) {
  txt <- list(
    paste0("    Ploidy level of ", x$name.p1, ":"),
    paste0("    Ploidy level of ", x$name.p2, ":"),
    paste0("    No. individuals:"),
    paste0("    No. markers:"),
    paste0("    Percentage of missing:"),
    paste0("    Percentage of redundant:")
  )
  n <- sapply(txt, nchar)
  for (i in 1:length(txt)) {
    txt[[i]] <- paste(txt[[i]], paste0(rep(" ", max(n) - n[i]), collapse = ""))
  }
  id <- is.na(x$geno.dose)
  cat("\n", txt[[1]], x$ploidy.p1)
  cat("\n", txt[[2]], x$ploidy.p2)
  cat("\n", txt[[3]], x$n.ind)
  cat("\n", txt[[4]], x$n.mrk)
  cat("\n ", txt[[5]], " (",   round(100*sum(id)/length(id),1), "%)", sep = "")
  if(!is.null(x$redundant))
    cat("\n ", txt[[6]], " (",   round(100*nrow(x$redundant)/sum(x$n.mrk, nrow(x$redundant)),1), "%)", sep = "")
  w <- table(x$chrom)
  w <- w[order(as.integer(gsub("[^0-9]", "", names(w))))]
  if (all(is.null(x$chrom)) || all(is.na(x$chrom)))
    cat("\n     No. markers per sequence: not available")
  else {
    cat("\n     ----------\n     No. markers per sequence:\n")
    print(data.frame(chrom = paste0("       ", names(w)), No.mrk = as.numeric(w)), row.names = FALSE)
  }
  cat("     ----------\n     No. of markers per dosage in both parents:\n")
  freq <- table(paste(x$dosage.p1,
                      x$dosage.p2, sep = "-"))
  d.temp <- matrix(unlist(strsplit(names(freq), "-")), ncol = 2, byrow = TRUE)
  d.temp <- data.frame(paste0("    ", d.temp[, 1]),
                       d.temp[, 2],
                       as.numeric(freq))
  colnames(d.temp) <- c(x$name.p1, x$name.p2, "freq")
  print(d.temp, row.names = FALSE)
}

#' @rdname read_geno_csv
#' @export
#' @importFrom graphics barplot layout mtext image legend
#' @importFrom grDevices colorRampPalette
plot.mappoly2.data <- function(x, thresh.line = NULL, ...)
{
  oldpar <- par(mar = c(5,4,1,2))
  on.exit(par(oldpar))
  if(is.null(thresh.line))
    thresh.line <- 0.05/length(x$mrk.names)
  freq <- table(paste(x$dosage.p1, x$dosage.p2, sep = "-"))
  d.temp <- matrix(unlist(strsplit(names(freq), "-")), ncol = 2, byrow = TRUE)
  type <- apply(d.temp, 1, function(x,ploidy.p1, ploidy.p2) paste0(sort(abs(abs(as.numeric(x)-(ploidy.p1/2))-(ploidy.p2/2))), collapse = ""),
                ploidy.p1 = x$ploidy.p1, ploidy.p2 = x$ploidy.p2)
  type.names <- names(table(type))
  mrk.dist <- as.numeric(freq)
  names(mrk.dist) <- apply(d.temp, 1 , paste, collapse = "-")
  layout(matrix(c(1,1,1,2,3,3,6,4,5), 3, 3), widths = c(1.2,3,.5), heights = c(1.5,2,3))
  barplot(mrk.dist, las = 2, #col = pal[match(type, type.names)],
          xlab = "Number of markers",
          ylab = "Dosage combination", horiz = TRUE)
  pval <- x$chisq.pval
  if(is.null(x$chisq.pval))
  {
    plot(0, 0, axes = FALSE, xlab = "", ylab = "", type = "n")
    text(x = 0, y = 0, labels = "No segregation test", cex = 2)
  } else{
    par(mar = c(1,1,1,2))
    par(xaxs = "i")
    plot(log10(pval), axes = FALSE, xlab = "", ylab = "", pch = 16,
         col = rgb(red = 0.25, green = 0.64, blue = 0.86, alpha = 0.3))
    axis(4, line = 1)
    mtext(text = bquote(log[10](P)), side = 4, line = 4, cex = .7)
    lines(x = c(0, length(x$mrk.names)), y = rep(log10(thresh.line),2), col = 2, lty = 2)
  }
  par(mar = c(5,1,0,2))
  pal <- c("black", colorRampPalette(c("#D73027", "#F46D43", "#FDAE61", "#FEE090",
                                       "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1",
                                       "#4575B4"))(x$ploidy.p1/2 + x$ploidy.p2/2 + 1))
  names(pal) <- c(-1:(x$ploidy.p1/2 + x$ploidy.p2/2))
  M <- as.matrix(x$geno.dose[x$mrk.names,])
  M[is.na(M)] <- -1
  image(x = 1:nrow(M), z = M, axes = FALSE, xlab = "",
        col = pal[as.character(sort(unique(as.vector(M))))], useRaster = TRUE)
  mtext(text = "Markers", side = 1, line = .4)
  mtext(text = "Individuals", side = 2, line = .2)
  par(mar = c(0,0,0,0))
  plot(0:10,0:10, type = "n", axes = FALSE, xlab = "", ylab = "")
  legend(0,10,
         horiz = FALSE,
         legend = c("missing", 0:(x$ploidy.p1/2 + x$ploidy.p2/2)),
         pch = 22,
         pt.cex = 3,
         pt.bg = pal, pt.lwd = 0,
         bty = "n", xpd = TRUE)
  if(!is.null(x$redundant)){
    par(mar = c(5,0,2,2))
    red = round(100*nrow(x$redundant)/(length(x$mrk.names)),1)
    mat = matrix(c(100-red, red), ncol = 1)
    w = barplot(mat, main = "",
                xlab = "", col = c(blues9[3],blues9[6]),
                axes = F, width = .5, border = NA, xlim = c(0,1))

    text(w, c((100-red)/2,   100 - red/2),  c(paste0(100 - red, " %"), paste0(red, " %")))
    mtext(text = "Unique vs. Redundant", line = -1, side = 4, cex = .8)
  }
  par(mfrow = c(1,1))
}
