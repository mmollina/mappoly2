#' Order Genetic Sequence in a Mapping Data Object
#'
#' This function orders the genetic sequence within a mapping data object. It can
#' operate on different types of genetic data (e.g., MDS or genomic) and allows
#' for customization of the ordering process through various parameters.
#'
#' @param x A mapping data object that contains genetic mapping information.
#' @param lg An optional vector specifying the linkage groups to be ordered.
#'           If NULL, all linkage groups in the data object are considered.
#' @param type The type of genetic data to be ordered, either 'mds' for
#'             Multi-Dimensional Scaling or 'genome' for genomic data.
#' @param p An optional parameter for the ordering function, specific to the
#'          method being used (e.g., MDS).
#' @param n An optional parameter for the ordering function, specific to the
#'          method being used.
#' @param ndim The number of dimensions to be used in the ordering process.
#'             This parameter is primarily relevant for MDS.
#' @param weight.exponent The exponent for weighting in the ordering process.
#' @param verbose A logical value indicating whether to print detailed output
#'                during the execution of the function.
#'
#' @return The function returns the modified mapping data object with updated
#'         order of the genetic sequence.
#'
#' @details The function iterates over the specified linkage groups (or all groups
#'          if none are specified) and applies either MDS or genomic ordering
#'          based on the 'type' parameter. Additional parameters like 'p', 'n',
#'          and 'ndim' are used to fine-tune the ordering process.
#'
#' @export
order_sequence <- function(x,
                           lg = NULL,
                           type = c("mds", "genome"),
                           p = NULL,
                           n = NULL,
                           ndim = 2,
                           weight.exponent = 2,
                           verbose = TRUE){
  y <- parse_lg_and_type(x,lg,type)
  for(i in y$lg){
    if(verbose) cat("  -->", i)
    if(y$type == "mds")
      x$maps[[i]][[y$type]]$order <- mds(x, x$maps[[i]][[y$type]]$mkr.names,p,n,ndim,weight.exponent,verbose)
    else
      x$maps[[i]][[y$type]]$order <- genome_order(x, x$maps[[i]][[y$type]]$mkr.names, verbose = TRUE)
    if(verbose) cat("\n")
  }
  return(x)
}

#' Estimates Loci Position Using Multidimensional Scaling
#'
#' This function estimates loci positions using the Multidimensional Scaling (MDS) method.
#' The method is adapted from the package \code{MDSmap} and is based on the approach
#' proposed by Preedy and Hackett (2016).
#'
#' @param x An object of class \code{mappoly2.sequence}.
#' @param mrk.id A vector of marker IDs to be used in the analysis.
#' @param p Integer; the smoothing parameter for the principal curve.
#'          If \code{NULL} (default), this will be determined using leave-one-out cross-validation.
#' @param n Vector of integers or strings containing loci to be omitted from the analysis.
#' @param ndim Integer; the number of dimensions to be considered in the MDS procedure (default = 2).
#' @param weight.exponent Integer; the exponent used in the LOD score values to weight the MDS procedure (default = 2).
#' @param verbose Boolean; if \code{TRUE} (default), displays information about the analysis.
#'
#' @return A list containing various components related to MDS results, including the input
#'         distance map, unconstrained MDS results, principal curve results, matrix of pairwise
#'         distances, data frame of loci positions, total length of the segment, vector of
#'         removed loci, scaling factor, and data frames for interpreting MDS plots.
#'
#' @details The function employs MDS to estimate the positions of loci on a genetic map.
#'          It includes options for excluding certain loci, adjusting the number of dimensions,
#'          and setting weighting factors for LOD scores. The results include detailed information
#'          about the estimated positions and relationships between loci.
#'
#' @references
#' Preedy, K. F., & Hackett, C. A. (2016). A rapid marker ordering approach for
#' high-density genetic linkage maps in experimental autotetraploid populations
#' using multidimensional scaling. _Theoretical and Applied Genetics_, 129(11),
#' 2117-2132. \doi{10.1007/s00122-016-2761-8}
#'
#' @importFrom smacof smacofSym
#' @importFrom princurve principal.curve
#' @importFrom stats runif
#' @author Marcelo Mollinari, adapted from MDSmap codes by Katharine F. Preedy
#' @export
mds <- function(x,
                mrk.id,
                p = NULL,
                n = NULL,
                ndim = 2,
                weight.exponent = 2,
                verbose = TRUE)
{
  rf.mat <- x$data$pairwise.rf$rec.mat[mrk.id, mrk.id]
  lod.mat <- x$data$pairwise.rf$lod.mat[mrk.id, mrk.id]
  o <- is.na(rf.mat)
  rf.mat[o] <- 1e-07
  lod.mat[o] <- 1e-07
  if(weight.exponent != 1)
    lod.mat <- lod.mat^weight.exponent
  diag(lod.mat) <- diag(rf.mat) <- NA
  locinames <- rownames(rf.mat)
  lodrf <- list(rf = rf.mat, lod = lod.mat, nloci = ncol(rf.mat), locinames = locinames)
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
  ord.obj <- list(smacofsym = smacofsym,
                  pc = pc1,
                  distmap = distmap,
                  lodmap = lodmap,
                  locimap = locimap,
                  length = max(estpos),
                  removed = n,
                  locikey = locikey,
                  meannnfit = nnfit$meanfit,
                  ndim = ndim)
  return(ord.obj)
}


#' Plot Multi-Dimensional Scaling (MDS) Maps
#'
#' This function visualizes the results of Multi-Dimensional Scaling (MDS) analysis
#' on genetic mapping data. It provides options to display either two-dimensional
#' or three-dimensional MDS plots based on the number of dimensions in the MDS object.
#'
#' @param x A genetic mapping data object containing MDS information.
#' @param lg A numeric value specifying the linkage group to be visualized.
#'           Defaults to the first linkage group.
#' @param D1lim Optional range for the first dimension of the MDS plot.
#' @param D2lim Optional range for the second dimension of the MDS plot.
#' @param D3lim Optional range for the third dimension of the MDS plot, applicable
#'              only for three-dimensional plots.
#' @param displaytext A logical value indicating whether to display text labels
#'                    on the MDS plot. Defaults to `FALSE`.
#'
#' @return The function does not return a value but generates a plot of the MDS
#'         analysis results for the specified linkage group.
#'
#' @details Depending on the dimensionality of the MDS object (`ndim`), this
#'          function calls either `plot_pcmap` for two-dimensional data or
#'          `plot_pcmap3d` for three-dimensional data. It visualizes the spatial
#'          arrangement of markers in the specified dimensions, providing insights
#'          into the genetic structure.
#'
#' @importFrom graphics plot text lines
#' @importFrom assertthat assert_that
#' @export
plot_mds <- function(x,
                     lg = 1,
                     D1lim = NULL,
                     D2lim = NULL,
                     D3lim = NULL,
                     displaytext = FALSE){
  assert_that(is.mds.ordered(x, lg))
  obj <- x$maps[[lg]][["mds"]]$order
  if(obj$ndim > 2){
    plot_pcmap3d(obj, D1lim, D2lim, D3lim, displaytext)
  } else {
    plot_pcmap(obj, D1lim, D2lim, displaytext)
  }
}


#' @author Katharine F. Preedy, \email{katharine.preedy@bioss.ac.uk}
#' @keywords internal
plot_pcmap <- function (x, D1lim = NULL, D2lim = NULL, displaytext = FALSE, ...)
{
  oldpar <- par(mfrow = c(1, 2))
  on.exit(par(oldpar))
  with(x, {
    if (displaytext  ==  TRUE) {
      labels = x$locikey$locus
    }
    else {
      labels = x$locikey$confplotno
    }
    graphics::plot(x$smacofsym$conf, type = "n", main = "MDS with principal curve",
                   xlim = D1lim, ylim = D2lim, xlab = "Dim 1", ylab = "Dim 2")
    text(x$smacofsym$conf, labels = labels, cex = 0.8)
    lines(pc)
    if (displaytext  ==  TRUE) {
      labels1 = x$locimap$locus
    }
    else {
      labels1 = x$locimap$confplotno
    }
    graphics::plot(x$locimap$position, x$locimap$nnfit, type = "n",
                   xlab = "Position", ylab = "nnfit", main = "nearest neighbour fits")
    text(x$locimap$position, x$locimap$nnfit, labels1)
  })
}


#' @author Katharine F. Preedy, \email{katharine.preedy@bioss.ac.uk}
#' @keywords internal
plot_pcmap3d <- function(x, D1lim = NULL, D2lim = NULL, D3lim = NULL, displaytext = FALSE, ...)
{
  oldpar <- par(mfrow = c(2, 2))
  on.exit(par(oldpar))
  with(x, {
    if (displaytext  ==  TRUE) {
      labels = x$locikey$locus
    }
    else {
      labels = x$locikey$confplotno
    }
    graphics::par(mfrow = c(2, 2))
    graphics::plot(x$smacofsym$conf[, "D1"], x$smacofsym$conf[,
                                                          "D2"], type = "n", main = "MDS with principal curve",
                   xlab = "Dimension 1", ylab = "Dimension 2", xlim = D1lim,
                   ylim = D2lim)
    text(x$smacofsym$conf[, "D1"], x$smacofsym$conf[, "D2"],
         labels = labels, cex = 0.8)
    lines(x$pc$s[, "D1"][x$pc$ord], x$pc$s[, "D2"][pc$ord])
    graphics::plot(x$smacofsym$conf[, "D1"], x$smacofsym$conf[,
                                                          "D3"], type = "n", main = "MDS with principal curve",
                   xlab = "Dimension 1", ylab = "Dimension 3", xlim = D1lim,
                   ylim = D3lim)
    text(x$smacofsym$conf[, "D1"], x$smacofsym$conf[, "D3"],
         labels = labels, cex = 0.8)
    lines(x$pc$s[, "D1"][pc$ord], x$pc$s[, "D3"][x$pc$ord])
    graphics::plot(x$smacofsym$conf[, "D2"], x$smacofsym$conf[,
                                                          "D3"], type = "n", main = "MDS with principal curve",
                   xlab = "Dimension 2", ylab = "Dimension 3", xlim = D2lim,
                   ylim = D3lim)
    text(x$smacofsym$conf[, "D2"], x$smacofsym$conf[, "D3"],
         labels = labels, cex = 0.8)
    lines(x$pc$s[, "D2"][x$pc$ord], pc$s[, "D3"][x$pc$ord])
    if (displaytext  ==  TRUE) {
      labels1 = x$locimap$locus
    }
    else {
      labels1 = x$locimap$confplotno
    }
    graphics::plot(x$locimap$position, x$locimap$nnfit, type = "n",
                   xlab = "Position", ylab = "nnfit", main = "nearest neighbour fits")
    text(x$locimap$position, x$locimap$nnfit, labels1)
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


#' Get the Genomic Position of Markers in a Sequence
#'
#' This internal function retrieves and orders the genomic position of markers in a given genetic sequence. It is primarily used within a larger analytical context.
#'
#' @param x An object of class \code{mappoly2.sequence}.
#' @param mrk.names A vector of marker names for which the genomic position is required.
#' @param verbose Logical; if TRUE, progress messages will be printed.
#'
#' @return Returns a data frame with the genomic position of the specified markers,
#'         ordered by chromosome and sequence position.
#'
#' @details The function checks for the availability of genomic position information
#'          in the provided data object. If available, it orders the markers based
#'          on their chromosome and sequence position. The function handles cases
#'          where only chromosome information or only sequence position information
#'          is available.
#'
#' @importFrom utils head
#' @importFrom assertthat assert_that
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
genome_order <- function(x, mrk.names, verbose = TRUE){
  assert_that(is.mappoly2.sequence(x))
  genome.pos <- x$data$genome.pos[mrk.names]
  chrom <- x$data$chrom[mrk.names]
  if(all(is.na(genome.pos))){
    if(all(is.na(chrom)))
      stop("No sequence or sequence position information found.")
    else{
      if (verbose) message("Ordering markers based on chromosome information")
      M <- data.frame(seq = chrom, row.names = mrk.names)
      M.out <- M[order(embedded_to_numeric(M[,1])),]
    }
  } else if(all(is.na(chrom))){
    if(all(is.na(genome.pos)))
      stop("No sequence or sequence position information found.")
    else{
      if (verbose) message("Ordering markers based on sequence position information")
      M <- data.frame(seq.pos = genome.pos, row.names = mrk.names)
      M.out <- M[order(embedded_to_numeric(M[,1])),]
    }
  } else{
    M <- data.frame(seq = chrom,
                    seq.pos = genome.pos,
                    row.names = mrk.names)
    M.out <- M[order(embedded_to_numeric(M[, "seq"]),
                     embedded_to_numeric(M[, "seq.pos"])),]
  }
  return(M.out)
}

