#' Create a sequence of markers
#'
#' Makes a sequence of markers based on an object of another class.
#'
#' @param x an object of one of the following classes:
#'     \code{mappoly.data}, \code{mappoly.map}, \code{mappoly.sequence},
#'     \code{mappoly.group},
#'     \code{mappoly.pcmap}, \code{mappoly.pcmap3d}, or \code{mappoly.geno.ord}
#'
#' @param arg one of the following input types:
#' \enumerate{
#' \item A string 'all': Generates a sequence with all markers from the raw data.
#' \item A string or a vector of strings 'chrX': Specifies a chromosome (with 'X'
#'       being the chromosome number). For unassigned markers, use 'chr0'.
#' \item A vector of integers: Indicates the position of markers in the original
#'       data set to be included in the sequence.
#' \item A vector of strings: Indicates the names or identifiers of the genetic
#'       markers
#' \item An integer: Represents a linkage group when xect is of class
#'       \code{mappoly.group}.
#' \item NULL: Applicable when \code{xect} belongs to one of the following
#'       classes: \code{mappoly.pcmap}, \code{mappoly.pcmap3d},
#'       \code{mappoly.unique.seq}, or \code{mappoly.geno.ord}.
#' }
#'
#' @param info.parent one of the following options:
#' \enumerate{
#' \item \code{'all'}{select all dosage combinations in both parents (default)}
#' \item \code{'P1'}{select informative markers parent 1}
#' \item \code{'P2'}{select informative markers parent 2}
#' }
#'
#' @param genomic.info An optional argument applicable only to \code{mappoly.group}
#' objects, which can be either NULL or a numeric combination of sequences
#' from genomic information used to create the sequences:
#' \enumerate{
#' \item NULL (default): Returns a sequence containing all markers as defined by
#' the grouping function.
#' \item 1: Returns a sequence with markers that match the intersection between
#' the grouping function and genomic information, considering the genomic information
#' sequence with the maximum number of matching markers for the group.
#' \item c(1, 2): Returns a sequence with markers that match the intersection between
#' the grouping function and genomic information, considering the two genomic
#' information sequences with the maximum number of matching markers for the group,
#' and so on.
#' }
#'
#' @param x an object of the class \code{mappoly.sequence}
#'
#' @param thresh.line position of a threshold line for p values of the segregation test (default = \code{0.05/n.mrk})
#'
#' @param ... currently ignored
#'
#' @return An object of class \code{mappoly2.sequence}, which is a
#'     list containing the following components:
#'     \item{mrk.names}{}
#'     \item{phases}{}
#'     \item{redundant}{}
#'     \item{data}{}
#'
#' @author Marcelo Mollinari (\email{mmollin@ncsu.edu}) and
#'         Gabriel Gesteira, (\email{gdesiqu@ncsu.edu})
#' @export
#' @importFrom assertthat assert_that
set_working_sequence <- function(x,
                                 lg = NULL,
                                 ch = NULL,
                                 mrk.names = NULL,
                                 ind.names = NULL,
                                 seq.names = NULL){
  assert_that(inherits(x, "screened"))
  ## If there is marker information, initiate a sequence with it
  if(!is.null(mrk.names)){
    return(initiate_working_sequence(x, mrk.names, ind.names, seq.names))
  }
  if(is.null(lg) & is.null(ch)){
    if(inherits(x, "grouped")){
      mrk.names <- split(names(x$linkage.groups$groups.snp), x$linkage.groups$groups.snp)
      seq.names <- paste0("Lg_", unique(x$linkage.groups$groups.snp))
      return(initiate_working_sequence(x, mrk.names, ind.names, seq.names))
    }
  } else if(!is.null(lg) & !is.null(ch)){
    assert_that(is.list(lg), msg = "'lg' should be a list")
    assert_that(is.list(ch), msg = "'ch' should be a list")
    assert_that(length(lg) == length(ch))
    mrk.names <- vector("list", length(lg))
    for(i in 1:length(lg)){
      mrk.names[[i]] <- mappoly2:::get_markers_from_grouped_and_chromosome(x, lg[[i]], ch[[i]])
    }
    return(initiate_working_sequence(x, mrk.names, ind.names, seq.names))
  } else if(!is.null(lg) & is.null(ch)){
    assert_that(is.list(lg), msg = "'lg' should be a list")
    mrk.names <- vector("list", length(lg))
    for(i in 1:length(lg)){
      mrk.names[[i]] <- mappoly2:::get_markers_from_grouped_sequence(x, lg[[i]])
    }
    return(initiate_working_sequence(x, mrk.names, ind.names, seq.names))
  } else if(is.null(lg) & !is.null(ch)){
    mrk.names <- vector("list", length(ch))
    assert_that(is.list(ch), msg = "'ch' should be a list")
    for(i in 1:length(ch)){
      mrk.names[[i]] <- rownames(mappoly2:::get_markers_from_chromosome(x, ch[[i]]))
      if(!is.null(x$initial.screened.rf))
        mrk.names[[i]] <- intersect(mrk.names[[i]], x$initial.screened.rf$mrk.names)
      else if(!is.null(x$screened.data$mrk.names)){
        mrk.names[[i]] <- intersect(mrk.names[[i]], x$screened.data$mrk.names)
      }
    }
    return(initiate_working_sequence(x, mrk.names, ind.names, seq.names))
  }
}

#' @importFrom graphics barplot layout mtext image legend
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices blues9
plot_sequence <- function(x, thresh.line = NULL, ...){
  oldpar <- par(mar = c(5,4,1,2))
  on.exit(par(oldpar))
  if(is.null(thresh.line))
    thresh.line <- 0.05/length(x$mrk.names)
  freq <- table(paste(x$data$dosage.p1[x$mrk.names], x$data$dosage.p2[x$mrk.names], sep = "-"))
  d.temp <- matrix(unlist(strsplit(names(freq), "-")), ncol = 2, byrow = TRUE)
  type <- apply(d.temp, 1, function(x,ploidy.p1, ploidy.p2) paste0(sort(abs(abs(as.numeric(x)-(ploidy.p1/2))-(ploidy.p2/2))), collapse = ""),
                ploidy.p1 = x$data$ploidy.p1, ploidy.p2 = x$data$ploidy.p2)
  type.names <- names(table(type))
  mrk.dist <- as.numeric(freq)
  names(mrk.dist) <- apply(d.temp, 1 , paste, collapse = "-")
  layout(matrix(c(1,1,1,2,3,3,6,4,5), 3, 3), widths = c(1.2,3,.5), heights = c(1.5,2,3))
  barplot(mrk.dist, las = 2, #col = pal[match(type, type.names)],
          xlab = "Number of markers",
          ylab = "Dosage combination", horiz = TRUE)
  pval <- x$data$chisq.pval[x$mrk.names]
  if(is.null(x$data$chisq.pval))
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
                                       "#4575B4"))(x$data$ploidy.p1/2 + x$data$ploidy.p2/2 + 1))
  names(pal) <- c(-1:(x$data$ploidy.p1/2 + x$data$ploidy.p2/2))
  M <- as.matrix(x$data$geno.dose[x$mrk.names,])
  M[is.na(M)] <- -1
  image(x = 1:nrow(M), z = M, axes = FALSE, xlab = "",
        col = pal[as.character(sort(unique(as.vector(M))))], useRaster = TRUE)
  mtext(text = "Markers", side = 1, line = .4)
  mtext(text = "Individuals", side = 2, line = .2)
  par(mar = c(0,0,0,0))
  plot(0:10,0:10, type = "n", axes = FALSE, xlab = "", ylab = "")
  legend(0,10,
         horiz = FALSE,
         legend = c("missing", 0:(x$data$ploidy.p1/2 + x$data$ploidy.p2/2)),
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

get_markers_from_chromosome <- function(x, arg){
  assert_that(has.chromosome.info(x))
  pattern <- "(ch|chr|CH|Chr|CHR|chrom|Chrom|Chromsome)"
  if(is.numeric(arg))
    ch.n.arg <- arg
  else{
    assert_that(all(is.character(arg)))
    assert_that(sum(grepl(pattern, arg, ignore.case = TRUE))  ==  length(arg))
    ch.n.arg <- mappoly2:::embedded_to_numeric(arg)
  }
  ch.n.dat <- mappoly2:::embedded_to_numeric(x$data$chrom)
  ch.id <- ch.n.dat%in%ch.n.arg
  mrk.names <- x$data$mrk.names[ch.id]
  chrom <- x$data$chrom[mrk.names]
  data.frame(chrom)
}

get_markers_from_grouped_sequence <- function(x, arg){
  inherits(x, "grouped")
  assert_that(is.numeric(arg))
  names(x$linkage.groups$groups.snp)[x$linkage.groups$groups.snp  %in%  arg]
}

get_markers_from_grouped_and_chromosome <- function(x, lg = NULL, ch = NULL){
  inherits(x, "grouped")
  assert_that(!is.null(lg) | !is.null(ch), msg = "Please provide a value for either 'lg' or 'ch'. Both cannot be left blank.")
  mrk.id.ch <- x$initial.sequence
  if(!is.null(ch))
    mrk.id.ch <- rownames(mappoly2:::get_markers_from_chromosome(x, ch))
  assert_that(length(mrk.id.ch) > 0)
  if(!is.null(lg))
    mrk.id.lg <- mappoly2:::get_markers_from_grouped_sequence(x, lg)
  return(intersect(mrk.id.ch, mrk.id.lg))
}

#' @export
set_initial_sequence <- function(x, arg){
  assert_that(inherits(x, "screened"))
  pattern <- "(ch|chr|CH|Chr|CHR|chrom|Chrom|Chromsome)"
  if(all(arg == "all"))
  {
    mrk.names <- x$screened.data$mrk.names

  } else if (all(is.character(arg)) &
             sum(grepl(pattern, arg, ignore.case = TRUE))  ==  length(arg) &
             all(!arg%in%rownames(x$data$geno.dose)))
  {
    if (all(is.na(x$data$chrom)))
      stop("There is no chromosome information.")
    ch.n.arg <- embedded_to_numeric(arg)
    ch.n.dat <- embedded_to_numeric(x$data$chrom)
    ch.id <- ch.n.dat%in%ch.n.arg
    mrk.names <- x$data$mrk.names[ch.id]
    chrom <- x$data$chrom[mrk.names]
    if (any(!is.na(x$data$genome.pos)))
      genome.pos <- x$data$genome.pos[mrk.names]
    ch_geno <- data.frame(chrom, genome.pos)
    sorted_ch_geno <- ch_geno[with(ch_geno, order(chrom, genome.pos)),]
    mrk.names <- rownames(sorted_ch_geno)
  } ## sequence with specific markers
  else if (all(is.character(arg)) & (length(arg)  ==  length(arg %in% x$data$mrk.names)))
  {
    mrk.names <- intersect(arg, x$data$mrk.names)
  }
  else if (is.vector(arg) && all(is.numeric(arg)))
  {
    assert_that(max(arg) <= x$data$n.mrk)
    mrk.names <- x$data$mrk.names[arg]
  }
  x$initial.sequence <- intersect(mrk.names, x$screened.data$mrk.names)
  class(x) <- unique(c(class(x), "initiated"))
  return(x)
}

initiate_working_sequence <- function(x,
                                      mrk.names,
                                      ind.names = NULL,
                                      seq.names = NULL){
  assert_that(inherits(x, "screened"))
  if(!is.list(mrk.names)){
    assert_that(is.character(mrk.names))
    mrk.names <- list(mrk.names)
  }
  if(is.null(ind.names)){
    ind.names <- vector("list", length(mrk.names))
    for(i in 1:length(ind.names)){
      ind.names[[i]] <- x$screened.data$ind.names
    }
  } else if(is.character(ind.names)){
    ind.names <- list(ind.names)
  }
  assert_that(length(ind.names) == length(mrk.names))
  x$working.sequences <- vector("list", length(mrk.names))
  if(is.null(seq.names))
    seq.names <- paste0("Seq_", 1:length(mrk.names))
  names(x$working.sequences) <- seq.names
  for(i in 1:length(mrk.names)){
    x$working.sequences[[i]] <- list(mrk.names = intersect(x$screened.data$mrk.names, mrk.names[[i]]),
                                     ind.names = intersect(x$screened.data$ind.names, ind.names[[i]]),
                                     order = list(mds = list(info = NULL,
                                                              phase = NULL),
                                                  genome = list(info = NULL,
                                                                 phase = NULL)))
  }
  class(x) <- unique(c(class(x), "ws"))
  return(x)
}

