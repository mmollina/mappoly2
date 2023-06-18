#' Create a sequence of markers
#'
#' Makes a sequence of markers based on an object of another class.
#'
#' @param input.obj an object of one of the following classes:
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
#' \item An integer: Represents a linkage group when input.object is of class
#'       \code{mappoly.group}.
#' \item NULL: Applicable when \code{input.object} belongs to one of the following
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
make_sequence <- function(input.obj,
                     arg = NULL,
                     info.parent = c("both", "p1", "p2"),
                     genomic.info = NULL,
                     phase = NULL) {
  info.parent <- match.arg(info.parent)
  if (is.mappoly2.data(input.obj))
  {
    pattern <- "(ch|chr|CH|Chr|CHR|chrom|Chrom|Chromsome)"
    ## Sequence with all markers
    if (all(arg  ==  "all"))
    {
      mrk.names <- input.obj$mrk.names
      out.dat <- input.obj
    } ## If chromosome informed
    else if (all(is.character(arg)) &
             sum(grepl(pattern, arg, ignore.case = TRUE))  ==  length(arg) &
             all(!arg%in%rownames(input.obj$geno.dose)))
    {
      if (all(is.na(input.obj$chrom)))
        stop("There is no chromosome information.")
      ch.n.arg <- embedded_to_numeric(arg)
      ch.n.dat <- embedded_to_numeric(input.obj$chrom)
      ch.id <- ch.n.dat%in%ch.n.arg
      mrk.names <- input.obj$mrk.names[ch.id]
      chrom <- input.obj$chrom[mrk.names]
      if (any(!is.na(input.obj$genome.pos)))
        genome.pos <- input.obj$genome.pos[mrk.names]
      ch_geno <- data.frame(chrom, genome.pos)
      sorted_ch_geno <- ch_geno[with(ch_geno, order(chrom, genome.pos)),]
      mrk.names <- rownames(sorted_ch_geno)
      out.dat <- subset_data(input.obj, select.mrk = mrk.names)
    } ## sequence with specific markers
    else if (all(is.character(arg)) & (length(arg)  ==  length(arg %in% input.obj$mrk.names)))
    {
      mrk.names <- intersect(arg, input.obj$mrk.names)
      out.dat <- subset_data(input.obj, select.mrk = mrk.names)
    }
    else if (is.vector(arg) && all(is.numeric(arg)))
    {
      assert_that(max(arg) <= input.obj$n.mrk)
      mrk.names <- input.obj$mrk.names[arg]
      out.dat <- subset_data(input.obj, select.mrk = mrk.names)
    }
    else stop("Invalid argument to select markers")
  }
  else if (is.mappoly2.sequence(input.obj))
  {
    return(make_sequence(input.obj$data, arg, info.parent))
  }
  else if (is.mappoly2.group(input.obj))
  {
    lgs.idx <- names(input.obj$groups.snp[input.obj$groups.snp  %in%  arg])
    if(is.null(genomic.info)){
      return(make_sequence(input.obj = input.obj$input.seq,
                      arg = lgs.idx))
    } else {
      assert_that(is.numeric(genomic.info))
      chrom <- input.obj$input.seq$data$chrom[lgs.idx]
      chrom.table <- sort(table(chrom, useNA = "always"), decreasing = TRUE)
      seq.group <- names(chrom)[chrom %in% names(chrom.table[genomic.info])]
      return(make_sequence(input.obj$input.seq$data, seq.group))
    }
  }
  else if (is.mappoly2.geno.ord(input.obj))
  {
    if(!is.null(arg))
      warning("Ignoring argument 'arg' and using the genome order instead.")
    return(make_sequence(input.obj$data, rownames(input.obj$ord)))
  }
  if (is.mappoly2.pcmap(input.obj) | is.mappoly2.pcmap3d(input.obj))
  {
    if(!is.null(arg))
      warning("Ignoring argument 'arg' and using the MDS order instead.")
    return(input.obj$mds.seq)
  }
  d.p1 <- input.obj$dosage.p1[mrk.names]
  d.p2 <- input.obj$dosage.p2[mrk.names]
  if(info.parent == "p1"){
    mrk.names <- mrk.names[d.p2 == 0 | d.p2 == input.obj$ploidy.p2]
    out.dat <- subset_data(input.obj, select.mrk = mrk.names)
  }
  else if(info.parent == "p2"){
    mrk.names <- mrk.names[d.p1 == 0 | d.p1 == input.obj$ploidy.p1]
    out.dat <- subset_data(input.obj, select.mrk = mrk.names)
  }
  structure(list(mrk.names = mrk.names,
                 phases = phase,
                 redundant = NULL,
                 data = out.dat),
            class = "mappoly2.sequence")
}

#' @rdname make_sequence
#' @export
print.mappoly2.sequence <- function(x, detailed = FALSE,  ...) {
  txt <- list(
    paste0("    Ploidy level of ", x$data$name.p1, ":"),
    paste0("    Ploidy level of ", x$data$name.p2, ":"),
    paste0("    No. individuals:"),
    paste0("    No. markers:"),
    paste0("    Percentage of missing:"),
    paste0("    Phases:"),
    paste0("      '--> Number of configurations:"),
    paste0("      '--> Percentage phased:"))
  n <- sapply(txt, nchar)
  for (i in 1:length(txt)) {
    txt[[i]] <- paste(txt[[i]], paste0(rep(" ", max(n) - n[i]), collapse = ""))
  }
  id <- is.na(x$data$geno.dose[x$mrk.names, ])
  cat("\n", txt[[1]], x$data$ploidy.p1)
  cat("\n", txt[[2]], x$data$ploidy.p2)
  cat("\n", txt[[3]], x$data$n.ind)
  cat("\n", txt[[4]], length(x$mrk.names))
  cat("\n ", txt[[5]], " ",   round(100*sum(id)/length(id),1), "%", sep = "")
  cat("\n", txt[[6]])
  if(is.null(x$phases)){
    cat("\n ", txt[[7]], " 0", sep = "")
    cat("\n ", txt[[8]], " 0%", sep = "")
  } else {
    cat("\n ", txt[[7]], " ", length(x$phases), sep = "")
    cat("\n ", txt[[8]], " ", nrow(x$phases[[1]]$p1), " (",   round(100*nrow(x$phases[[1]]$p1)/length(x$mrk.names),1), "%)", sep = "")
  }
  w <- table(x$data$chrom[x$mrk.names], useNA = "always")
  w <- w[order(as.integer(gsub("[^0-9]", "", names(w))))]
  names(w)[is.na(names(w))] <- "NoCrh"
  if (all(is.null(x$data$chrom[x$mrk.names])) || all(is.na(x$data$chrom[x$mrk.names])))
    cat("\n     No. markers per sequence: not available")
  else {
    cat("\n     ----------\n     No. markers per sequence:\n")
    print(data.frame(chrom = paste0("       ", names(w)), No.mrk = as.numeric(w)), row.names = FALSE)
  }

  if(detailed){
    cat("     ----------\n     No. of markers per dosage in both parents:\n")
    freq <- table(paste(x$data$dosage.p1[x$mrk.names],
                        x$data$dosage.p2[x$mrk.names], sep = "-"))
    d.temp <- matrix(unlist(strsplit(names(freq), "-")), ncol = 2, byrow = TRUE)
    d.temp <- data.frame(paste0("    ", d.temp[, 1]),
                         d.temp[, 2],
                         as.numeric(freq))
    colnames(d.temp) <- c(x$data$name.p1, x$data$name.p2, "freq")
    print(d.temp, row.names = FALSE)
  }
}

#' @rdname make_sequence
#' @export
#' @importFrom graphics barplot layout mtext image legend
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices blues9
plot.mappoly2.sequence <- function(x, thresh.line = NULL, ...)
{
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



