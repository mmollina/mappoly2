#' @export
#' @importFrom graphics barplot layout mtext image legend
#' @importFrom grDevices colorRampPalette
plot.mappoly2 <- function(x,
                          what = c("raw",
                                   "screened",
                                   "initiated",
                                   "pairwise",
                                   "group"),
                          type = c("rf", "lod"),
                          ord = NULL,
                          rem = NULL,
                          main.text = NULL,
                          index = FALSE,
                          fact = 1, ...)
{
  what <- match.arg(what)
  if (what == "raw"){
    assert_that(inherits(x, "data"))
    plot_data(x, text = "Raw data", col = "darkred")
  }
  else if(what == "screened"){
    assert_that(inherits(x, "screened"))
      plot_data(x,
                text = "All data - screened",
                col =  "#07607e",
                mrk.id = x$screened.data$mrk.names,
                ind.id = x$screened.data$ind.names)
  }
  else if(what == "initiated"){
    assert_that(inherits(x, "initiated"))
    plot_data(x,
              text = "Initial sequence - screened",
              col =  "#07607e",
              mrk.id = x$initial.sequence,
              ind.id = x$screened.data$ind.names)
  }
  else if(what == "pairwise"){
    assert_that(inherits(x, "pairwise"))
    plot_rf_matrix(x$pairwise,
                   type = type,
                   ord = ord,
                   rem = rem,
                   main.text = main.text,
                   index = index,
                   fact = fact)
  } else if(what == "group"){
    plot_group(x$linkage.groups)
  }
}


plot_data <- function(x, text, col, mrk.id = NULL, ind.id = NULL, ...){
  oldpar <- par(mar = c(5,4,1,2))
  on.exit(par(oldpar))
  if(is.null(mrk.id))
    mrk.id <- x$data$mrk.names
  if(is.null(ind.id))
    ind.id <- x$data$ind.names
  freq <- table(paste(x$data$dosage.p1[mrk.id],
                      x$data$dosage.p2[mrk.id],
                      sep = "-"))
  d.temp <- matrix(unlist(strsplit(names(freq), "-")), ncol = 2, byrow = TRUE)
  type <- apply(d.temp, 1, function(x,ploidy.p1, ploidy.p2) paste0(sort(abs(abs(as.numeric(x)-(ploidy.p1/2))-(ploidy.p2/2))), collapse = ""),
                ploidy.p1 = x$data$ploidy.p1, ploidy.p2 = x$data$ploidy.p2)
  type.names <- names(table(type))
  mrk.dist <- as.numeric(freq)
  names(mrk.dist) <- apply(d.temp, 1 , paste, collapse = "-")
  layout(matrix(c(1,1,1,2,3,3,6,4,5), 3, 3), widths = c(1.2,3,.5), heights = c(1.5,2,3))
  barplot(mrk.dist, las = 2,
          xlab = "Number of markers",
          ylab = "Dosage combination", horiz = TRUE)
  pval <- x$data$QAQC.values$markers[,"chisq.pval"]
  names(pval) <- rownames(x$data$QAQC.values$markers)
  if(is.null(pval))
  {
    plot(0, 0, axes = FALSE, xlab = "", ylab = "", type = "n")
    text(x = 0, y = 0, labels = "No segregation test", cex = 2)
  } else{
    pval <- pval[mrk.id]
    par(mar = c(1,1,1,2))
    par(xaxs = "i")
    plot(log10(pval), axes = FALSE, xlab = "", ylab = "", pch = 16,
         col = rgb(red = 0.25, green = 0.64, blue = 0.86, alpha = 0.3))
    axis(4, line = 1)
    mtext(text = bquote(log[10](P)), side = 4, line = 4, cex = .7)
  }
  par(mar = c(5,1,0,2))
  pal <- c("black", colorRampPalette(c("#D73027", "#F46D43", "#FDAE61", "#FEE090",
                                       "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1",
                                       "#4575B4"))(x$data$ploidy.p1/2 + x$data$ploidy.p2/2 + 1))
  names(pal) <- c(-1:(x$data$ploidy.p1/2 + x$data$ploidy.p2/2))
  M <- as.matrix(x$data$geno.dose[mrk.id,ind.id])
  M[is.na(M)] <- -1
  image(x = 1:nrow(M), z = M, axes = FALSE, xlab = "",
        col = pal[as.character(sort(unique(as.vector(M))))], useRaster = TRUE)
  mtext(text = "Markers", side = 1, line = .4)
  mtext(text = text, side = 1, line = 2, col = col, cex = 0.5)
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
  if(!is.null(x$data$redundant)){
    par(mar = c(5,0,2,2))
    red = round(100*nrow(x$data$redundant)/(length(x$data$mrk.names) + nrow(x$data$redundant)),1)
    if(length(red) == 0){
      plot(x = 0, y = 0, type = "n", axes = F, xlab = "", ylab = "")
      text(x = 0, y = 0, "No redundant markers", srt=90)
    }
    else {
      mat = matrix(c(100-red, red), ncol = 1)
      w = barplot(mat, main = "",
                  xlab = "", col = c(blues9[3],blues9[6]),
                  axes = F, width = .5, border = NA, xlim = c(0,1))

      text(w, c((100-red)/2,   100 - red/2),  c(paste0(100 - red, " %"), paste0(red, " %")))
      mtext(text = "Unique vs. Redundant", line = -1, side = 4, cex = .8)
    }
  }
  par(mfrow = c(1,1))
}
plot_rf_matrix <- function(x, type = c("rf", "lod"), ord = NULL, rem = NULL,
                           main.text = NULL, index = FALSE, fact = 1, ...){
  type <- match.arg(type)
  #if(ncol(x$rec.mat) > 1000)
  #  fact = 3
  #if(ncol(x$rec.mat) > 5000)
  #  fact = 5
  #if(ncol(x$rec.mat) > 10000)
  #  fact = 10
  if(is.mappoly2.sequence(ord))
    ord <- ord$mrk.names
  if(type  ==  "rf"){
    w <- x$rec.mat
    if(!is.null(ord))
    {
      w <- w[ord,ord]
    }
    if(!(is.null(rem) || sum(colnames(x$rec.mat)%in%rem)  ==  0))
    {
      o <- which(colnames(x$rec.mat)%in%rem)
      w <- w[-o,-o]
    }
    if(fact > 1)
      w <- aggregate_matrix(w, fact)
    if(is.null(main.text))
      main.text <- "Recombination fraction matrix"
    col.range  <-
      na.omit(rev(fields::tim.colors())[1:(ceiling(128 * max(x$rec.mat, na.rm = TRUE)) + 1)])
    brks <- NULL
  } else if(type  ==  "lod")
  {
    w <- x$lod.mat
    if(!is.null(ord))
    {
      w <- w[ord,ord]
    }
    if(!(is.null(rem) || sum(colnames(x$rec.mat)%in%rem)  ==  0))
    {
      o <- which(colnames(x$rec.mat)%in%rem)
      w <- w[-o,-o]
    }
    if(fact > 1)
      w <- aggregate_matrix(w, fact)
    w[w < 1e-4] <- 1e-4
    w <- log10(w)
    if(is.null(main.text))
      main.text <- "log(LOD) Score matrix"
    col.range <- na.omit(fields::tim.colors()[1:(ceiling(128 * max(x$lod.mat, na.rm = TRUE)) + 1)])
    col.range <- col.range[ceiling(seq(1, length(col.range), length.out = 10))]
    brks <- seq(min(w, na.rm = TRUE), max(w, na.rm = TRUE), length.out = 11)
    brks <- round(exp(brks/log10(exp(1))),1)
  } else stop("Invalid matrix type.")

  fields::image.plot(
    w,
    col = col.range,
    lab.breaks = brks,
    main = main.text,
    useRaster = FALSE,
    axes = FALSE
  )
  if(ncol(w) < 100)
    ft <- .7
  else
    ft <- 100/ncol(w)
  if(index)
    text(x = seq(0,1, length.out = ncol(w)), y = seq(0,1, length.out = ncol(w)),
         labels = colnames(w), cex = ft)
}
plot_group <- function(x, ...) {
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

