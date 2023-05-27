#' Estimates loci position using Multidimensional Scaling
#'
#' Estimates loci position using the Multidimensional Scaling method proposed by
#' \cite{Preedy and Hackett (2016)}. This code is adapted from
#' the package \code{MDSmap}, available under the GNU General Public License,
#' Version 3, at \url{https://CRAN.R-project.org/package=MDSMap}.
#'
#' @param input.mat an object of class \code{mappoly.input.matrix}.
#'
#' @param p integer. The smoothing parameter for the principal curve.
#'   If \code{NULL} (default), this will be determined using leave-one-out cross validation.
#'
#' @param n vector of integers or strings containing loci to be omitted from the analysis.
#'
#' @param ndim integer. The number of dimensions to be considered in the multidimensional scaling procedure (default = 2).
#'
#' @param weight.exponent integer. The exponent used in the LOD score values to weight the
#'        MDS procedure (default = 2).
#'
#' @param verbose boolean. If \code{TRUE} (default), displays information about the analysis.
#'
#' @param x an object of class \code{mappoly.mds}.
#'
#' @param ... currently ignored.
#'
#' @return A list containing:
#' \item{M}{the input distance map.}
#' \item{sm}{the unconstrained MDS results.}
#' \item{pc}{the principal curve results.}
#' \item{distmap}{a matrix of pairwise distances between
#' loci, where the columns are in the estimated order.}
#' \item{locimap}{a data frame of the loci containing the name
#' and position of each locus in order of increasing distance.}
#' \item{length}{integer representing the total length of the segment.}
#' \item{removed}{a vector of the names of loci removed from the analysis.}
#' \item{scale}{the scaling factor from the MDS.}
#' \item{locikey}{a data frame showing the number associated with each
#' locus name for interpreting the MDS configuration plot.}
#' \item{confplotno}{a data frame showing locus names associated
#' with each number on the MDS configuration plots.}
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu} mostly adapted from MDSmap
#'         codes, written by Katharine F. Preedy, \email{katharine.preedy@bioss.ac.uk}
#'
#' @references
#'  Preedy, K. F., & Hackett, C. A. (2016). A rapid marker ordering approach for
#'  high-density genetic linkage maps in experimental autotetraploid populations
#'  using multidimensional scaling. _Theoretical and Applied Genetics_, 129(11),
#'  2117-2132. \doi{10.1007/s00122-016-2761-8}
#'
#' @importFrom smacof smacofSym
#' @importFrom princurve principal.curve
#' @importFrom stats runif
#' @importFrom utils read.csv write.csv
#' @export mds
mds <- function(input.mat,
                      p = NULL,
                      n = NULL,
                      ndim = 2,
                      weight.exponent = 2,
                      verbose = TRUE)
{
  o <- is.na(input.mat$rec.mat)
  input.mat$rec.mat[o] <- 1e-07
  input.mat$lod.mat[o] <- 1e-07
  if(weight.exponent != 1)
    input.mat$lod.mat <- input.mat$lod.mat^weight.exponent
  diag(input.mat$lod.mat) <- diag(input.mat$rec.mat) <- NA
  locinames <- rownames(input.mat$rec.mat)
  lodrf <- list(rf = input.mat$rec.mat, lod = input.mat$lod.mat, nloci = ncol(input.mat$rec.mat), locinames = locinames)
  confplotno <- 1:lodrf$nloci
  if(!is.null(n)){
    if(!is.numeric(n))n <- which(lodrf$locinames%in%n)
    r <- lodrf$rf[-n,-n]
    lod <- lodrf$lod[-n,-n]
    confplotno <- confplotno[-n]
  } else {
    r <- lodrf$rf
    lod <- lodrf$lod
  }
  M <- imf_h(r)/100
  nloci = length(confplotno)
  smacofsym <- smacof::smacofSym(M,ndim = ndim,weightmat = lod,itmax = 100000)
  pc1 <- princurve::principal_curve(smacofsym$conf,maxit = 150,spar = p,smoother = "smooth_spline")
  scale <- sum(smacofsym$delta)/sum(smacofsym$dhat)
  # Configuration dissim are based on the normalized observed diss - dhat.
  # True observed dissimilarities are delta
  maporder <- pc1$ord
  estpos <- pc1$lambda[maporder]*scale*100
  # gives the estimated length from the beginning of the line
  rownames <- lodrf$locinames[maporder]
  distmap <- outer(maporder,maporder,Vectorize(function(i,j)M[i,j]))
  lodmap <- outer(maporder,maporder, Vectorize(function(i,j)lod[i,j]))
  rownames(distmap) <- rownames;colnames(distmap) <- rownames
  rownames(lodmap) <- rownames;colnames(lodmap) <- rownames
  if(!is.null(n))  {
    locikey <- data.frame(locus = lodrf$locinames[-n],confplotno = confplotno)
  } else {
    locikey <- data.frame(locus = lodrf$locinames,confplotno = confplotno)
  }
  nnfit <- calc.nnfit(distmap,lodmap,estpos)
  locimap <- data.frame(confplotno = confplotno[maporder],locus = locikey$locus[maporder],position = estpos,nnfit = nnfit$pointfits,row.names = 1:nloci)
  if(!is.null(n)) {
    removedloci <- data.frame(n,lodrf$locinames[n],row.names = NULL)
  } else {
    removedloci <- n
  }
  map <- list(smacofsym = smacofsym,pc = pc1,distmap = distmap,lodmap = lodmap,locimap = locimap,length = max(estpos),removed = n,locikey = locikey,meannnfit = nnfit$meanfit)
  if(verbose)
  {
    cat(paste('Stress:', round(map$smacofsym$stress,5)))
    cat(paste('\nMean Nearest Neighbour Fit:', round(map$meannnfit,5)))
  }
  map$mds.seq <- make_sequence(input.mat$input.seq, as.character(map$locimap$locus))
  if(ndim  ==  2) {
    return(structure(map, class = "mappoly2.pcmap"))
  } else {
    return(structure(map, class = "mappoly2.pcmap3d"))
  }
}

#' @rdname mds
#' @export
print.mappoly2.pcmap <- function(x, ...)
{
  cat("\nThis is an object of class 'mappoly.mds'")
  cat("\nNumber of markers: ", nrow(x$locimap))
  cat("\nNumber of dimensions used: ", x$sm$ndim)
  cat("\nStress: ", x$smacofsym$stress)
  cat("\nMean Nearest Neighbour Fit:", x$meannnfit)
}

#' @rdname mds
#' @export
print.mappoly2.pcmap3d <- function(x, ...)
{
  cat("\nThis is an object of class 'mappoly.mds'")
  cat("\nNumber of markers: ", nrow(x$locimap))
  cat("\nNumber of dimensions used: ", x$sm$ndim)
  cat("\nStress: ", x$smacofsym$stress)
  cat("\nMean Nearest Neighbour Fit:", x$meannnfit)
}

#' @author Katharine F. Preedy, \email{katharine.preedy@bioss.ac.uk}
#' @export
plot.mappoly2.pcmap <- function (x, D1lim = NULL, D2lim = NULL, displaytext = FALSE, ...)
{
  oldpar <- par(mfrow = c(1, 2))
  on.exit(par(oldpar))
  with(x, {
    if (displaytext  ==  TRUE) {
      labels = locikey$locus
    }
    else {
      labels = locikey$confplotno
    }
    graphics::plot(smacofsym$conf, type = "n", main = "MDS with principal curve",
                   xlim = D1lim, ylim = D2lim, xlab = "Dim 1", ylab = "Dim 2")
    text(smacofsym$conf, labels = labels, cex = 0.8)
    lines(pc)
    if (displaytext  ==  TRUE) {
      labels1 = locimap$locus
    }
    else {
      labels1 = locimap$confplotno
    }
    graphics::plot(locimap$position, locimap$nnfit, type = "n",
                   xlab = "Position", ylab = "nnfit", main = "nearest neighbour fits")
    text(locimap$position, locimap$nnfit, labels1)
  })
}

#' @author Katharine F. Preedy, \email{katharine.preedy@bioss.ac.uk}
#' @export
plot.mappoly2.pcmap3d <- function(x, D1lim = NULL, D2lim = NULL, D3lim = NULL, displaytext = FALSE, ...)
{
  oldpar <- par(mfrow = c(2, 2))
  on.exit(par(oldpar))
  with(x, {
    if (displaytext  ==  TRUE) {
      labels = locikey$locus
    }
    else {
      labels = locikey$confplotno
    }
    graphics::par(mfrow = c(2, 2))
    graphics::plot(smacofsym$conf[, "D1"], smacofsym$conf[,
                                                          "D2"], type = "n", main = "MDS with principal curve",
                   xlab = "Dimension 1", ylab = "Dimension 2", xlim = D1lim,
                   ylim = D2lim)
    text(smacofsym$conf[, "D1"], smacofsym$conf[, "D2"],
         labels = labels, cex = 0.8)
    lines(pc$s[, "D1"][pc$ord], pc$s[, "D2"][pc$ord])
    graphics::plot(smacofsym$conf[, "D1"], smacofsym$conf[,
                                                          "D3"], type = "n", main = "MDS with principal curve",
                   xlab = "Dimension 1", ylab = "Dimension 3", xlim = D1lim,
                   ylim = D3lim)
    text(smacofsym$conf[, "D1"], smacofsym$conf[, "D3"],
         labels = labels, cex = 0.8)
    lines(pc$s[, "D1"][pc$ord], pc$s[, "D3"][pc$ord])
    graphics::plot(smacofsym$conf[, "D2"], smacofsym$conf[,
                                                          "D3"], type = "n", main = "MDS with principal curve",
                   xlab = "Dimension 2", ylab = "Dimension 3", xlim = D2lim,
                   ylim = D3lim)
    text(smacofsym$conf[, "D2"], smacofsym$conf[, "D3"],
         labels = labels, cex = 0.8)
    lines(pc$s[, "D2"][pc$ord], pc$s[, "D3"][pc$ord])
    if (displaytext  ==  TRUE) {
      labels1 = locimap$locus
    }
    else {
      labels1 = locimap$confplotno
    }
    graphics::plot(locimap$position, locimap$nnfit, type = "n",
                   xlab = "Position", ylab = "nnfit", main = "nearest neighbour fits")
    text(locimap$position, locimap$nnfit, labels1)
  })
}

#'@author Katharine F. Preedy, \email{katharine.preedy@bioss.ac.uk}
#'@keywords internal
calc.nnfit <- function (distmap, lodmap, estmap){
  pointfits <- unlist(lapply(1:dim(distmap)[2], calc.nnfit.loci,
                             distmap = distmap, lodmap = lodmap, estmap = estmap))
  fit <- sum(pointfits)
  list(fit = fit, pointfits = pointfits, meanfit = mean(pointfits))
}

#'@author Katharine F. Preedy, \email{katharine.preedy@bioss.ac.uk}
#'@keywords internal
calc.nnfit.loci <- function (loci, distmap, lodmap, estmap){
  nns <- get.nearest.informative(loci, lodmap)
  obs <- distmap[loci, nns]
  est <- estmap[loci] - estmap[nns]
  nn.fit <- sum(abs(obs - est))
  nn.fit
}

#'@author Katharine F. Preedy, \email{katharine.preedy@bioss.ac.uk}
#'@keywords internal
get.nearest.informative <- function (loci, lodmap){
  neighbours <- NULL
  if (loci > 1) {
    locileft <- lodmap[loci, (loci - 1):1]
    if (length(which(locileft != 0)) > 0)
      neighbours <- loci - min(which(locileft != 0))
  }
  if (loci < dim(lodmap)[2]) {
    lociright <- lodmap[loci, (loci + 1):dim(lodmap)[2]]
    if (length(which(lociright != 0)) > 0)
      neighbours <- c(neighbours, loci + min(which(lociright != 0)))
  }
  neighbours
}


#' Get the genomic position of markers in a sequence
#'
#' This functions gets the genomic position of markers in a sequence and
#' return an ordered data frame with the name and position of each marker
#' @param input.seq a sequence object of class \code{mappoly.sequence}
#' @param verbose if \code{TRUE} (default), the current progress is shown; if
#'     \code{FALSE}, no output is produced
#' @param x 	an object of the class mappoly.geno.ord
#' @param ... 	currently ignored
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export get_genome_order
#' @importFrom utils head
get_genome_order <- function(input.seq, verbose = TRUE){
  assert_that(is.mappoly2.sequence(input.seq))
  genome.pos <- input.seq$data$genome.pos[input.seq$mrk.names]
  chrom <- input.seq$data$chrom[input.seq$mrk.names]
  if(all(is.na(genome.pos))){
    if(all(is.na(chrom)))
      stop("No sequence or sequence position information found.")
    else{
      if (verbose) message("Ordering markers based on chromosome information")
      M <- data.frame(seq = chrom, row.names = input.seq$mrk.names)
      M.out <- M[order(embedded_to_numeric(M[,1])),]
    }
  } else if(all(is.na(chrom))){
    if(all(is.na(genome.pos)))
      stop("No sequence or sequence position information found.")
    else{
      if (verbose) message("Ordering markers based on sequence position information")
      M <- data.frame(seq.pos = genome.pos, row.names = input.seq$mrk.names)
      M.out <- M[order(embedded_to_numeric(M[,1])),]
    }
  } else{
    M <- data.frame(seq = chrom,
                    seq.pos = genome.pos,
                    row.names = input.seq$mrk.names)
    M.out <- M[order(embedded_to_numeric(M[, "seq"]),
                     embedded_to_numeric(M[, "seq.pos"])),]
  }
  structure(list(data = input.seq$data, ord = M.out), class = "mappoly2.geno.ord")
}

#' @rdname get_genome_order
#' @export
print.mappoly2.geno.ord <- function(x, ...){
  print(head(x$ord))
}

#' @rdname get_genome_order
#' @export
#' @importFrom ggplot2 ggplot aes geom_point xlab ylab theme
plot.mappoly2.geno.ord <- function(x, ...){
  seq <- seq.pos <- NULL
  w <- x$ord
  w$seq <- as.factor(w$seq)
  w$seq = with(w, reorder(seq))
  ggplot2::ggplot(w,
                  ggplot2::aes(x = seq.pos, y = seq, group = as.factor(seq))) +
    ggplot2::geom_point(ggplot2::aes(color = as.factor(seq)), shape = 108, size = 5, show.legend = FALSE) +
    ggplot2::xlab("Position") +
    ggplot2::ylab("Chromosome")
}


