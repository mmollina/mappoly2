#' Assign Markers to Linkage Groups
#'
#' This function identifies linkage groups of markers by utilizing the results
#' from two-point (pairwise) analysis.
#'
#' @param input.mat A \code{mappoly2.rf.matrix} object containing the
#' recombination fraction and LOD score matrices.
#'
#' @param input.seq A \code{mappoly2.sequence} object containing the markers
#' to be grouped. If \code{NULL} (default), all markers in \code{input.mat} are used.
#'
#' @param expected.groups The number of expected linkage groups (e.g. chromosomes)
#' for the species, if known.
#'
#' @param inter A logical indicating whether to plot a dendrogram highlighting the
#' expected groups before continuing. Default is \code{TRUE}.
#'
#' @param comp.mat A logical indicating whether to display a comparison between the
#' reference-based and linkage-based groupings if chromosome information is available.
#' Default is \code{FALSE}.
#'
#' @param LODweight A logical indicating whether the clustering should be weighted by
#' the square of the LOD score. Default is \code{FALSE}.
#'
#' @return An object of class \code{mappoly2.group}, which is a list containing:
#' \itemize{
#' \item{hc.snp}{A list with information related to the UPGMA grouping method.}
#' \item{expected.groups}{The number of expected linkage groups.}
#' \item{groups.snp}{The assigned groups for each marker.}
#' \item{seq.vs.grouped.snp}{A comparison between the genomic group information
#' (when available) and the groups provided by \code{group}.}
#' }
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}#'
#' @importFrom graphics abline pie
#' @importFrom stats as.dendrogram as.dist cutree hclust lm predict quantile rect.hclust
#' @importFrom dendextend color_branches
#' @export group

group <- function(input.data = NULL,
                  use.genome.only = FALSE,
                  expected.groups = NULL,
                  inter = TRUE,
                  comp.mat = FALSE,
                  LODweight = FALSE)
{
  ## checking for correct object
  if(use.genome.only){


    .seq_skeleton()

  }

  assert_that(mappoly2:::has.mappoly2.rf(input.data), )
  MSNP <- input.seq$pairwise$rec.mat[input.seq$mrk.names, input.seq$mrk.names]
  mn <- input.seq$data$chrom[input.seq$mrk.names]
  mn[is.na(mn)] <- "NoChr"
  dimnames(MSNP) <- list(mn, mn)
  diag(MSNP) <- 0
  MSNP[is.na(MSNP)] <- .5
  if(LODweight){
    Mlod <- input.seq$pairwise$lod.mat^2
    Mlod[is.na(Mlod)] <- 10e-5
    hc.snp <- hclust(as.dist(MSNP), method="ward.D2", members=apply(Mlod, 1, mean, na.rm = TRUE))
  } else {
    hc.snp <- hclust(as.dist(MSNP), method = "average")
  }
  ANSWER <- "flag"
  if(interactive() && inter)
  {
    dend.snp <- as.dendrogram(hc.snp)
    while(substr(ANSWER, 1, 1) != "y" && substr(ANSWER, 1, 1) != "yes" && substr(ANSWER, 1, 1) != "Y" && ANSWER  != "")
    {
      if(is.null(expected.groups))
        expected.groups <- as.numeric(readline("Enter the number of expected groups: "))
      dend1 <- dendextend::color_branches(dend.snp, k = expected.groups)
      plot(dend1, leaflab = "none")
      z <- rect.hclust(hc.snp, k = expected.groups, border = "red")
      groups.snp  <- cutree(tree = hc.snp, k = expected.groups)
      xy <- sapply(z, length)
      xt <- as.numeric(cumsum(xy)-ceiling(xy/2))
      yt <- .1
      points(x = xt, y = rep(yt, length(xt)), cex = 6, pch = 20, col = "lightgray")
      text(x = xt, y = yt, labels = pmatch(xy, table(groups.snp, useNA = "ifany")), adj = .5)
      ANSWER <- readline("Enter 'Y/n' to proceed or update the number of expected groups: ")
      if(substr(ANSWER, 1, 1)  ==  "n" | substr(ANSWER, 1, 1)  ==  "no" | substr(ANSWER, 1, 1)  ==  "N")
        stop("Function halted.")
      if(substr(ANSWER, 1, 1) != "y" && substr(ANSWER, 1, 1) != "yes" && substr(ANSWER, 1, 1) != "Y" && ANSWER  != "")
        expected.groups <- as.numeric(ANSWER)
    }
  }
  if(is.null(expected.groups))
    stop("Inform the 'expected.groups' or use 'inter = TRUE'")
  # Distribution of SNPs into linkage groups
  seq.vs.grouped.snp <- NULL
  if(all(unique(mn)  ==  "NoChr") & comp.mat)
  {
    comp.mat <- NA
    seq.vs.grouped.snp <- NA
    warning("There is no physical reference to generate a comparison matrix")
  }
  groups.snp  <- cutree(tree = hc.snp, k = expected.groups)
  if(comp.mat){
    seq.vs.grouped.snp <- matrix(0, expected.groups, length(na.omit(unique(mn))),
                                 dimnames = list(1:expected.groups, na.omit(unique(mn))))
    for(i in 1:expected.groups)
    {
      x <- table(names(which(groups.snp == i)))
      seq.vs.grouped.snp[i,names(x)] <- x
    }
    idtemp2 <- unique(apply(seq.vs.grouped.snp, 1, which.max))
    idtemp2 <- c(idtemp2, setdiff(1:(ncol(seq.vs.grouped.snp)-1), idtemp2))
    seq.vs.grouped.snp <- cbind(seq.vs.grouped.snp[,idtemp2])
    cnm <- colnames(seq.vs.grouped.snp)
    colnames(seq.vs.grouped.snp) <- cnm
  } else {
    seq.vs.grouped.snp <- NULL
  }
  names(groups.snp) <- input.seq$mrk.names
  input.seq$linkage.groups <- list(hc.snp = hc.snp,
                                   expected.groups = expected.groups,
                                   groups.snp = groups.snp,
                                   seq.vs.grouped.snp = seq.vs.grouped.snp)
  return(input.seq)
}

#' @export
print_mappoly2_group <- function(x, detailed = TRUE, ...) {
  ## criteria
  cat("\n       - Number of linkage groups:  ", length(unique(x$groups.snp)), "\n")
  cat("       - Markers per linkage groups: \n")
  w <- table(x$groups.snp, useNA = "ifany")
  w <- data.frame(group = names(w), n_mrk = as.numeric(w), row.names = NULL)
  mappoly2:::print_matrix(mat = w, 8, row.names = FALSE)
  cat("\n")
  ## printing summary
  if(!is.null(x$seq.vs.grouped.snp)){
    mappoly2:::print_matrix(mat = x$seq.vs.grouped.snp, 8)
  }
}

#' @export
plot_mappoly2_group <- function(x, ...) {
  dend <- as.dendrogram(x$hc.snp)
  dend1 <- dendextend::color_branches(dend, k = x$expected.groups)
  plot(dend1, leaflab = "none")
  z <- rect.hclust(x$hc.snp, k = x$expected.groups, border = "red")
  xy <- sapply(z, length)
  xt <- as.numeric(cumsum(xy)-ceiling(xy/2))
  yt <- .1
  points(x = xt, y = rep(yt, length(xt)), cex = 6, pch = 20, col = "lightgray")
  text(x = xt, y = yt, labels = pmatch(xy, table(x$groups.snp, useNA = "ifany")), adj = .5)
}
