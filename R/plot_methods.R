#' Plot Data from mappoly2.data Object
#'
#' Visualizes data from a `mappoly2.data` object, which may also have "screened" and "pairwise.rf" classes. This function provides various plotting options including recombination frequency matrices, screened data, marker density, and raw data, depending on the type of data available in the object and the selected options.
#'
#' @param x A `mappoly2.data` object containing genetic mapping data. The object can also be of class "screened" and "pairwise.rf".
#' @param type A character string specifying the type of plot to generate. It can be one of "rf" (recombination frequency), "screened", "density", or "raw". Default is "rf".
#' @param chrom Optional; a vector of chromosome numbers or names to subset the data before plotting. If `NULL`, all chromosomes are considered.
#' @param ... Additional arguments, not used in this method
#' @details
#' The function can visualize data in four different ways based on the `type` parameter:
#' - "rf": Plots a recombination frequency matrix, available for objects with the "pairwise.rf" class.
#' - "screened": Shows screened data for objects with the "screened" class.
#' - "density": Generates a density plot of markers, requiring genome position information.
#' - "raw": Displays raw data from the `mappoly2.data` object.
#' The function handles data selection based on the provided chromosome information (`chrom`).
#'
#' @return The function does not return a value but generates a plot.
#'
#' @examples
#' plot(B2721)
#'
#' @importFrom graphics barplot layout mtext image legend
#' @importFrom grDevices colorRampPalette
#' @importFrom CMplot CMplot
#' @export
plot.mappoly2.data<-function(x,
                             type = c("rf", "screened", "density", "raw"),
                             chrom = NULL,
                             ...)
{
  mrk.id <- NULL
  if(!is.null(chrom)){
    mrk.id <- x$mrk.names[get_mrk_indices_from_chrom(x, chrom)]
  }
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  type <- match.arg(type)
  if(has.mappoly2.rf(x) & type == "rf"){
    assert_that(has.mappoly2.rf(x))
    if(!is.null(mrk.id))
      mrk.id <- Reduce(intersect, list(x$screened.data$mrk.names,
                                       mrk.id, colnames(x$pairwise.rf$rec.mat)))
    else
      mrk.id <- Reduce(intersect, list(x$screened.data$mrk.names,
                                       colnames(x$pairwise.rf$rec.mat)))
    plot_rf_matrix_one(x$pairwise.rf,
                       fact = ceiling(ncol(x$pairwise.rf$rec.mat)/1000),
                       ord = mrk.id)
  }
  else if (type == "density"){
    assert_that(has.mappoly2.screened(x))
    assert_that(data.has.genome.info(x))
    if(!is.null(mrk.id))
      mrk.id <- intersect(x$screened.data$mrk.names, mrk.id)
    else
      mrk.id <- x$screened.data$mrk.names
    u <- data.frame(SNP = mrk.id,
                    Chromosome = embedded_to_numeric(x$chrom[mrk.id]),
                    Position = x$genome.pos[mrk.id])
    CMplot::CMplot(u,type = "p", plot.type = "d", file.output=FALSE)
  }
  else if(((!has.mappoly2.rf(x) & has.mappoly2.screened(x)) | type == "screened") & type != "raw"){
    assert_that(has.mappoly2.screened(x))
    if(!is.null(mrk.id))
      mrk.id <- intersect(x$screened.data$mrk.names, mrk.id)
    else
      mrk.id <- x$screened.data$mrk.names
    plot_data(x, text = "Screened data", col = "darkblue",
              mrk.id = mrk.id,
              ind.id = x$screened.data$ind.names)
  }
  else
    plot_data(x, text = "Raw data", col = "darkred", mrk.id = mrk.id)
}

plot_data <- function(x,
                      text,
                      col,
                      mrk.id = NULL,
                      ind.id = NULL,
                      ...){
  oldpar <- par(mar = c(5,4,1,2))
  on.exit(par(oldpar))
  if(is.null(mrk.id))
    mrk.id <- x$mrk.names
  if(is.null(ind.id))
    ind.id <- x$ind.names
  freq <- table(paste(x$dosage.p1[mrk.id],
                      x$dosage.p2[mrk.id],
                      sep = "-"))
  d.temp <- matrix(unlist(strsplit(names(freq), "-")), ncol = 2, byrow = TRUE)
  type <- apply(d.temp, 1, function(x,ploidy.p1, ploidy.p2) paste0(sort(abs(abs(as.numeric(x)-(ploidy.p1/2))-(ploidy.p2/2))), collapse = ""),
                ploidy.p1 = x$ploidy.p1, ploidy.p2 = x$ploidy.p2)
  type.names <- names(table(type))
  mrk.dist <- as.numeric(freq)
  names(mrk.dist) <- apply(d.temp, 1 , paste, collapse = "-")
  layout(matrix(c(1,1,1,2,3,3,6,4,5), 3, 3), widths = c(1.2,3,.5), heights = c(1.5,2,3))
  barplot(mrk.dist, las = 2,
          xlab = "Number of markers",
          ylab = "Dosage combination", horiz = TRUE)
  pval <- x$QAQC.values$markers[,"chisq.pval"]
  names(pval) <- rownames(x$QAQC.values$markers)
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
                                       "#4575B4"))(x$ploidy.p1/2 + x$ploidy.p2/2 + 1))
  names(pal) <- c(-1:(x$ploidy.p1/2 + x$ploidy.p2/2))
  M <- as.matrix(x$geno.dose[mrk.id,ind.id])
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
         legend = c("missing", 0:(x$ploidy.p1/2 + x$ploidy.p2/2)),
         pch = 22,
         pt.cex = 3,
         pt.bg = pal, pt.lwd = 0,
         bty = "n", xpd = TRUE)
  if(!is.null(x$redundant)){
    par(mar = c(5,0,2,2))
    red = round(100*nrow(x$redundant)/(length(x$mrk.names) + nrow(x$redundant)),1)
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



#' Plot Recombination Fraction Matrices for Genetic Maps
#'
#' This function visualizes the recombination fraction matrices for genetic maps.
#' It allows the user to plot these matrices for different types of genetic data
#' (e.g., MDS, genome, custom) and for specified linkage groups.
#'
#' @param x A \code{mappoly2.data} object that contains recombination fraction information.
#' @param lg Optional vector specifying the linkage groups for which the recombination
#'           fraction matrices should be plotted. If NULL, matrices for all linkage
#'           groups are plotted.
#' @param type The type of genetic data to be visualized. Can be 'mds', 'genome',
#'             or 'custom'. This parameter determines how the matrices are processed
#'             and displayed.
#' @param fact A numeric factor used for scaling or aggregating the matrix data.
#'             Defaults to 1 (no scaling).
#'
#' @return The function does not return a value but generates a series of plots,
#'         each representing the recombination fraction matrix for a specified linkage
#'         group.
#'
#' @details The function processes the genetic mapping data based on the specified
#'          'type' and 'lg' parameters. It then uses `plot_rf_matrix_one` to plot
#'          individual matrices for each linkage group. The visualization helps in
#'          understanding the recombination patterns and genetic distances between
#'          markers.
#'
#' @importFrom graphics par
#' @importFrom assertthat assert_that
#' @export
plot_rf_matrix <- function(x,
                           lg = NULL,
                           type = c("mds", "genome", "custom"),
                           fact = 1){
  y <- parse_lg_and_type(x,lg,type)
  mrk.id <- get_markers_from_ordered_sequence(x, y$lg, y$type)
  op <- par(mfrow = optimal_layout(length(y$lg)), pty = "s")
  on.exit(par(op))
  for(i in 1:length(mrk.id)){
    plot_rf_matrix_one(x$data$pairwise.rf,
                       ord = mrk.id[[i]],
                       main.text = paste(names(x$maps[i]), y$type, sep = "-"),
                       fact = fact)
  }
}

plot_rf_matrix_one <- function(x,
                               type = c("rf", "lod"),
                               ord = NULL,
                               rem = NULL,
                               main.text = NULL,
                               index = FALSE,
                               fact = 1, ...){
  type <- match.arg(type)
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

#' @export
plot.mappoly2.group <- function(x, ...) {
  op <- par(lwd=1.5)
  on.exit(par(op))
  dend <- as.dendrogram(x$hc.snp)
  dend1 <- dendextend::color_branches(dend, k = x$expected.groups, col = mp_pal(x$expected.groups))
  plot(dend1, leaflab = "none")
  z <- list(names(x$groups.snp))
  if(x$expected.groups != 1){
    par(lwd=4)
    z <- rect.hclust(x$hc.snp, k = x$expected.groups, border = "black")
    par(lwd = 1)
  }
  xy <- sapply(z, length)
  xt <- as.numeric(cumsum(xy)-ceiling(xy/2))
  yt <- .1
  points(x = xt, y = rep(yt, length(xt)), cex = 6, pch = 20, col = "lightgray")
  text(x = xt, y = yt, labels = pmatch(xy, table(x$groups.snp, useNA = "ifany")), adj = .5)
}


#' prepare maps for plot
#' @param void internal function to be documented
#' @keywords internal
prepare_map <- function(x,
                        ploidy.p1, ploidy.p2,
                        name.p1, name.p2,
                        dosage.p1, dosage.p2,
                        alt=NULL, ref=NULL){
  ## Gathering marker positions
  map <- cumsum(imf_h(c(0, x$hmm.phase[[1]]$rf)))
  names(map) <- rownames(x$hmm.phase[[1]]$p1)
  ## Gathering phases
  ph.p1 <- x$hmm.phase[[1]]$p1
  ph.p2 <- x$hmm.phase[[1]]$p2
  colnames(ph.p1) <- paste0("p1.", 1:ploidy.p1)
  colnames(ph.p2) <- paste0("p2.", 1:ploidy.p2)
  if(is.null(ref))
  {
    ph.p1[ph.p1 == 1] <- ph.p2[ph.p2 == 1] <- "A"
    ph.p1[ph.p1 == 0] <- ph.p2[ph.p2 == 0] <- "B"
  } else {
    for(i in names(map)){
      ph.p1[i, ph.p1[i,] == 1] <- alt[i]
      ph.p1[i, ph.p1[i,] == 0] <- ref[i]
      ph.p2[i, ph.p2[i,] == 1] <- alt[i]
      ph.p2[i, ph.p2[i,] == 0] <- ref[i]
    }
  }
  d.p1 <- dosage.p1[names(map)]
  d.p2 <- dosage.p2[names(map)]
  list(ploidy.p1 = ploidy.p1,
       ploidy.p2 = ploidy.p2,
       name.p1 = name.p1,
       name.p2 = name.p2,
       map = map,
       ph.p1 = ph.p1,
       ph.p2 = ph.p2,
       d.p1 = d.p1,
       d.p2 = d.p2)
}


#' Plot Genetic Map
#'
#' This function visualizes a genetic map for a specified linkage group. It supports various types of genetic maps and offers customization options for the display.
#'
#' @param x An object representing genetic mapping data.
#' @param lg The linkage group to be visualized, default is the first linkage group (lg = 1).
#' @param type The type of genetic map to process, either "mds" or "genome".
#' @param parent Specifies which parent's data to use in the visualization.
#'               Options are "p1p2" (both parents), "p1" (first parent), or "p2" (second parent).
#' @param left.lim The left limit for the plotting area, default is 0.
#' @param right.lim The right limit for the plotting area, default is Inf.
#' @param phase Logical; if TRUE, phases are included in the plot.
#' @param mrk.names Logical; if TRUE, marker names are displayed on the plot.
#' @param plot.dose Logical; if TRUE, doses are plotted.
#' @param homolog.names.adj Adjustment for homolog names in the plot.
#' @param cex Character expansion size for text in the plot.
#' @param xlim The limits for the x-axis. Can be set to NULL for automatic adjustment.
#' @param main The main title for the plot.
#' @param ... Additional graphical parameters.
#'
#' @return The function does not return a value but generates a plot of the genetic map.
#'
#' @details The function creates a detailed plot of a genetic map for a given linkage group.
#'          It can display various features such as phases, marker names, and doses, and allows
#'          for customization of the plot's appearance.
#'
#' @importFrom grDevices rgb blues9
#' @importFrom graphics rect
#' @export
plot_map <- function(x, lg = 1, type = c("mds", "genome"),
                     parent = c("p1p2", "p1", "p2"),
                     left.lim = 0, right.lim = Inf,
                     phase = TRUE, mrk.names = FALSE,
                     plot.dose = TRUE, homolog.names.adj = 3,
                     cex = 1, xlim = NULL, main = "",...) {
  type <- match.arg(type)
  y <- parse_lg_and_type(x,lg,type)
  parent <- match.arg(parent)
  assert_that(is.mapped.sequence(x, y$lg, y$type, parent),
              msg = "Requested map is not estimated")
  assert_that(length(y$lg) ==1 & is.numeric(lg))

  v <- detect_hmm_est_map(x)
  u <- apply(v[parent,,,drop=FALSE],1,all)
  h <- names(u)[1:2][!u[1:2]]
  if(length(h) == 1)
    assert_that(u[type], msg = paste(h, "order has not been computed for", parent))
  else
    assert_that(u[type], msg = paste(h[1], "and", h[2],"orders have not been computed for", parent))


  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  map.info <- prepare_map(x$maps[[y$lg]][[y$type]][[parent]],
                          x$data$ploidy.p1, x$data$ploidy.p2,
                          x$data$name.p1, x$data$name.p2,
                          x$data$dosage.p1, x$data$dosage.p2,
                          x$data$alt, x$data$ref)
  if(any(map.info$ph.p1 == "B")){
    var.col <- c(A = "black", B = "darkgray")
  } else {
    var.col <- c(A = "#008000", T = "#FF0000", C = "#0000FF", G = "#FFFF00")
  }
  ploidy <- max(c(map.info$ploidy.p1, map.info$ploidy.p2))
  x <- map.info$map
  lab <- names(x)
  zy <- seq(0, 0.6, by = 0.12)
  zy.p1 <- zy[1:map.info$ploidy.p1] +1.8 + (0.3 * ((map.info$ploidy.p2/2)-1))
  zy.p2 <- zy[1:map.info$ploidy.p2] + 1.1
  pp <- map.info$ph.p1
  pq <- map.info$ph.p2
  d.p1 <- map.info$d.p1
  d.p2 <- map.info$d.p2
  x1 <- abs(left.lim - x)
  x2 <- abs(right.lim - x)
  id.left <- which(x1 == min(x1))[1]
  id.right <- rev(which(x2 == min(x2)))[1]
  par(mai = c(1,0.15,0,0), mar = c(4.5,homolog.names.adj,1,2))
  curx <- x[id.left:id.right]
  layout(mat  = matrix(c(2,4,1,3), ncol = 2), heights = c(10, 1), widths = c(1, 10))
  #layout(mat  = matrix(c(4,2,3, 1), ncol = 2), heights = c(2, 10), widths = c(1, 10))
  if(is.null(xlim)){
    xlim <- range(curx)
  }
  max.y <- 4.0
  plot(x = curx,
       y = rep(.5,length(curx)),
       type = "n" ,
       ylim = c(.25, max.y),
       axes = FALSE,
       xlab = "Distance (cM)",
       ylab = "",
       xlim = xlim)
  lines(c(x[id.left], x[id.right]), c(.5, .5), lwd = 15, col = "gray")
  points(x = curx,
         y = rep(.5,length(curx)),
         xlab = "", ylab = "",
         pch = "|", cex = 1.5,
         ylim = c(0,2))
  axis(side = 1)
  ##Parent1
  x1 <- seq(x[id.left], x[id.right], length.out = length(curx))
  x.control <- diff(x1[1:2])/2
  if(length(x1) < 150)
    x.control <- x.control * .8
  if(length(x1) < 100)
    x.control <- x.control * .8
  if(length(x1) < 75)
    x.control <- x.control * .8
  if(length(x1) < 50)
    x.control <- x.control * .8
  if(length(x1) < 25)
    x.control <- x.control * .8
  for(i in 1:map.info$ploidy.p2)
  {
    lines(range(x1), c(zy.p2[i], zy.p2[i]), lwd = 12, col = "gray")
    y1 <- rep(zy.p2[i], length(curx))
    pal <- var.col[pq[id.left:id.right,i]]
    rect(xleft = x1 - x.control,
         ybottom = y1 -.035,
         xright = x1 + x.control,
         ytop = y1 +.035,
         col = pal,
         border = NA)
  }
  #connecting allelic variants to markers
  for(i in 1:length(x1))
    lines(c(curx[i], x1[i]), c(0.575, zy.p2[1]-.05), lwd = 0.2)
  ####
  if(plot.dose){
    y <- zy.p2[map.info$ploidy.p2]+0.1 - ((map.info$ploidy.p2/2 - 1)*0.005)+d.p2[id.left:id.right]/20
    y.l <- zy.p2[map.info$ploidy.p2]+0.1 - ((map.info$ploidy.p2/2 - 1)*0.005)+ c(0:map.info$ploidy.p2)[id.left:id.right]/20
    #text(x = min(x1) - 1, y = mean(y.l), "doses", srt = 90)
    for(i in 1:length(y.l)){
      text(x = min(x1)-1,y=y.l[i],i-1, cex = .7)
      lines(range(x1), c(y.l[i], y.l[i]), lwd = .5, col = "gray")
    }

    points(x = x1,
           y = y,
           col = "darkgray",
           #col = d.col[as.character(d.p2[id.left:id.right])],
           pch = 19, cex = .7)
  }
  ##Parent2
  for(i in 1:map.info$ploidy.p1)
  {
    lines(range(x1), c(zy.p1[i], zy.p1[i]), lwd = 12, col = "gray")
    y1 <- rep(zy.p1[i], length(curx))
    pal <- var.col[pp[id.left:id.right,i]]
    rect(xleft = x1 - x.control,
         ybottom = y1 -.035,
         xright = x1 + x.control,
         ytop = y1 +.035,
         col = pal,
         border = NA)
  }
  ####
  if(plot.dose){
    y <- zy.p1[map.info$ploidy.p1]+0.1 -((map.info$ploidy.p1/2 - 1)*0.005)+d.p1[id.left:id.right]/20
    y.l <- zy.p1[map.info$ploidy.p1]+0.1 -((map.info$ploidy.p1/2 - 1)*0.005)+c(0:map.info$ploidy.p1)[id.left:id.right]/20
    #text(x = min(x1) - 1, y = mean(y.l), "doses", srt = 90)
    for(i in 1:length(y.l)){
      text(x = min(x1)-1,y=y.l[i],i-1, cex = .7)
      lines(range(x1), c(y.l[i], y.l[i]), lwd = .5, col = "gray")
    }
    points(x = x1,
           y = y,
           col = "darkgray",
           #col = d.col[as.character(d.p1[id.left:id.right])],
           pch = 19, cex = .7)
  }
  if(mrk.names)
    text(x = x1,
         y = rep(max(y)+.1, length(x1)),
         labels = names(curx),
         srt = 90, adj = 0, cex = cex *.6)
  par(mar = c(4.5,1,1,0), xpd = TRUE)
  plot(x = 0,
       y = 0,
       type = "n" ,
       axes = FALSE,
       ylab = "",
       xlab = "",
       ylim = c(.25, max.y))

  mtext(text = main, side = 2, at = mean(c(zy.p2, zy.p2)), line = -1, font = 4, cex = cex , adj = c(0,0))
  mtext(text = map.info$name.p2, side = 4, at = mean(zy.p2), line = -1, font = 4)
  for(i in 1:map.info$ploidy.p2)
    mtext(colnames(map.info$ph.p2)[i], line = 1, at = zy.p2[i], side = 4, las = 2, cex = 0.7 * cex)
  mtext(text = map.info$name.p1, side = 4, at = mean(zy.p1), line = -1, font = 4)
  for(i in 1:map.info$ploidy.p1)
    mtext(colnames(map.info$ph.p1)[i],  line = 1, at = zy.p1[i], side = 4, las = 2, cex = 0.7 * cex)
  par(mar = c(0,0,0,0), xpd = FALSE)
  plot(x = curx,
       y = rep(.5,length(curx)),
       type = "n" ,
       axes = FALSE,
       xlab = "",
       ylab = "")
  if(any(map.info$ph.p1 == "B")){
    legend("topleft", legend = c("A", "B"),
           fill  = c(var.col), #title = "Variants",
           box.lty = 0, bg = "transparent", ncol = 6)
  } else {
    legend("topleft", legend = c("A", "T", "C", "G", "-"),
           fill  = c(var.col, "white"),# title = "Nucleotides",
           box.lty = 0, bg = "transparent", ncol = 6)
  }
}

#' Plot Physical vs. Genetic Distance
#'
#' This function creates scatterplots to compare physical distance (in Mbp) with genetic
#' distance (in cM). It accepts either a single object or a list of objects of class
#' \code{mappoly.map}, plotting the relationship for each map provided.
#'
#' @param x A single object or a list of objects of class \code{mappoly.map}.
#' @param type Character vector indicating the type of genetic map to be used for the
#'             analysis. Options include "mds", "genome", or "custom". Default is c("mds", "genome").
#' @param parent Character vector specifying the parent or parents to be considered in the
#'               analysis. Options are "p1p2" (both parents), "p1" (first parent), or "p2"
#'               (second parent). Default is c("p1p2", "p1", "p2").
#' @param same.ch.lg Logical; if \code{TRUE}, only scatterplots for chromosomes and linkage
#'                   groups with the same number are displayed. Default is \code{FALSE}.
#' @param alpha Numeric; transparency factor for SNP points in the scatterplot. Default is 1/5.
#' @param size Numeric; size of the SNP points in the scatterplot. Default is 3.
#'
#' @return A ggplot object representing the scatterplot(s) of physical vs. genetic distances.
#'
#' @details The function generates scatterplots to visually compare the physical distance
#'          (measured in megabase pairs, Mbp) and the genetic distance (measured in centiMorgans, cM)
#'          for each map. This helps in understanding the relationship between physical and genetic
#'          distances in genetic studies.
#'
#' @importFrom ggplot2 ggplot geom_point facet_wrap facet_grid labs theme_bw theme element_text
#' @importFrom assertthat assert_that
#' @export
plot_genome_vs_map <- function(x,
                               type = c("mds", "genome"),
                               parent = c("p1p2", "p1", "p2"),
                               same.ch.lg = FALSE,
                               alpha = 1/2,
                               size = 2){
  assert_that(is.mappoly2.sequence(x))
  type <- match.arg(type)
  parent <- match.arg(parent)
  has.map <- logical(length(x$maps))
  names(has.map) <- names(x$maps)
  for(i in names(x$maps)){
    has.map[i] <- is.mapped.sequence(x, i, type, parent)
  }
  if(!all(has.map)){
    stop("maps are no available for \n-->  type:", type, "\n-->  parent:", parent)
  }
  w <- lapply(x$maps[has.map], function(y) y[[type]])
  geno.vs.map <- NULL
  for(i in 1:length(w)){
    LG <- genomic.pos <- map.pos <- NULL
    mrk.names <- rownames(w[[i]][[parent]]$hmm.phase[[1]]$p1)
    geno.vs.map <- rbind(geno.vs.map,
                         data.frame(mrk.names = mrk.names,
                                    map.pos = cumsum(imf_h(c(0, w[[i]][[parent]]$hmm.phase[[1]]$rf))),
                                    genomic.pos = x$data$genome.pos[mrk.names]/1e6,
                                    LG = as.factor(i),
                                    chr = as.factor(x$data$chrom[mrk.names])))
  }
  geno.vs.map$chr <- factor(geno.vs.map$chr, levels = sort(levels(geno.vs.map$chr)))
  num_colors <- length(unique(geno.vs.map$LG))  # Number of unique levels in LG
  color_palette <- mp_pal(num_colors)  # Get the right number of colors from your palette
  if(same.ch.lg){
    p <- ggplot2::ggplot(geno.vs.map, ggplot2::aes(genomic.pos, map.pos)) +
      ggplot2::geom_point(alpha = alpha, ggplot2::aes(colour = LG), size = size) +
      ggplot2::scale_colour_manual(values = color_palette) +  # Apply custom color palette
      ggplot2::facet_wrap(~LG, nrow = floor(sqrt(length(w)))) +
      ggplot2::labs(subtitle = "Linkage group", x = "Genome position (Mbp)", y = "Map position (cM)") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     legend.position = "none", plot.subtitle = ggplot2::element_text(hjust = 0.5))
  } else {
    p <- ggplot2::ggplot(geno.vs.map, ggplot2::aes(genomic.pos, map.pos)) +
      ggplot2::geom_point(alpha = alpha, ggplot2::aes(colour = LG), size = size) +
      ggplot2::scale_colour_manual(values = color_palette) +  # Apply custom color palette
      ggplot2::facet_grid(LG~chr) +
      ggplot2::labs(x = "Genome position (Mbp)", y = "Map position (cM)") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "none")
  }

  # Print the plot
  print(p)

}

#' Plot a List of Genetic Maps
#'
#' This function visualizes a list of genetic maps from a mappoly2 sequence object. It supports both horizontal and vertical orientations and allows for customization of plot colors.
#'
#' @param x A mappoly2 sequence object containing genetic map data.
#' @param horiz Logical; if TRUE, the maps are plotted horizontally, otherwise vertically.
#' @param type The type of genetic maps to be plotted, either "mds" or "genome".
#' @param parent Specifies which parent's data to use in the visualization.
#'               Options are "p1p2" (both parents), "p1" (first parent), or "p2" (second parent).
#' @param col The color used for plotting the maps, default is "lightgray".
#'            Can be a vector of colors to apply different colors to each map.
#'
#' @return The function does not return a value but generates a plot or a series of plots
#'         representing the genetic maps. It invisibly returns a data frame containing
#'         marker positions and linkage group information.
#'
#' @details The function iterates over the linkage groups in the provided mappoly2 sequence
#'          object and creates a plot (or a series of plots) showing the genetic map(s)
#'          for the specified linkage groups. The plot orientation can be set to either
#'          horizontal or vertical, and the color of the plots can be customized.
#'
#' @importFrom graphics plot axis
#' @export
plot_map_list <- function(x, horiz = TRUE,
                          type = c("mds", "genome"),
                          parent = c("p1p2", "p1", "p2"),
                          col = "lightblue"){
  assert_that(is.mappoly2.sequence(x))
  type <- match.arg(type)
  parent <- match.arg(parent)
  has.map <- logical(length(x$maps))
  names(has.map) <- names(x$maps)
  for(i in names(x$maps)){
    has.map[i] <- is.mapped.sequence(x, i, type, parent)
  }
  if(!all(has.map)){
    stop("maps are no available for \n-->  type:", type, "\n-->  parent:", parent)
  }
  w <- lapply(x$maps[has.map], function(y) y[[type]])
  if(length(col) == 1)
    col <- rep(col, length(w))
  z <- NULL
  max.dist <- max(sapply(w, function(x) sum(imf_h(x[[parent]]$hmm.phase[[1]]$rf))))
  if(horiz){
    plot(0,
         xlim = c(0, max.dist),
         ylim = c(0,length(w)+1),
         type = "n", axes = FALSE,
         xlab = "Map position (cM)",
         ylab = "Linkage groups")
    axis(1)
    for(i in 1:length(w)){
      d <- cumsum(imf_h(c(0, w[[i]][[parent]]$hmm.phase[[1]]$rf)))
      z <- rbind(z, data.frame(mrk = rownames(w[[i]][[parent]]$hmm.phase[[1]]$p1),
                               LG = names(w)[i], pos = d))
      plot_one_map(d, i = i, horiz = TRUE, col = col[i])
    }
    axis(2, at = 1:length(w), labels = names(w), lwd = 0, las = 2)
  } else{
    plot(0,
         ylim = c(-max.dist, 0),
         xlim = c(0,length(w)+1),
         type = "n", axes = FALSE,
         ylab = "Map position (cM)",
         xlab = "Linkage groups")
    x <- axis(2, labels = FALSE, lwd = 0)
    axis(2, at = x, labels = abs(x))
    for(i in 1:length(w)){
      d <- cumsum(imf_h(c(0, w[[i]][[parent]]$hmm.phase[[1]]$rf)))
      z <- rbind(z, data.frame(mrk = rownames(w[[i]][[parent]]$hmm.phase[[1]]$p1),
                               LG = names(w)[i], pos = d))
      plot_one_map(d, i = i, horiz = FALSE, col = col[i])
    }
    axis(3, at = 1:length(w), labels = names(w), lwd = 0, las = 2)
  }
  invisible(z)
}


plot_one_map<-function(x,
                       i = 0,
                       horiz = FALSE,
                       col = "lightgray")
{
  if(horiz)
  {
    rect(xleft = x[1], ybottom = i-0.25,
         xright = tail(x,1), ytop = i+0.25,
         col = col)
    for(j in 1:length(x))
      lines(x = c(x[j], x[j]), y = c(i-0.25, i+0.25), lwd = .5)
  } else {
    x <- -rev(x)
    rect(xleft = i-0.25, ybottom = x[1],
         xright = i+0.25, ytop = tail(x,1),
         col = col)
    for(j in 1:length(x))
      lines(y = c(x[j], x[j]), x = c(i-0.25, i+0.25), lwd = .5)
  }
}

#' Plot MDS vs. Genome
#'
#' This function generates a plot comparing MDS (Multi-Dimensional Scaling) positions to genome positions for genetic markers in a `mappoly2.sequence` object. It uses functions from the `mappoly2` package to extract marker data and `ggplot2` for visualization.
#'
#' @param x A `mappoly2.sequence` object containing marker and map data. This object typically results from mapping or marker analysis processes in the `mappoly2` package.
#' @param alpha The transparency level of the points in the plot, with 1 being fully opaque and 0 being fully transparent. Defaults to 1/2.
#' @param size The size of the points in the plot. Defaults to 2.
#'
#' @return A ggplot object representing the MDS versus genome position plot. Each linkage group is displayed in a separate panel.
#'
#' @importFrom ggplot2 ggplot geom_point facet_wrap labs theme_bw theme element_text
#'
#' @export
plot_mds_vs_genome <- function(x,
                               alpha = 1/2,
                               size = 2){
  d <- lg <- y <- NULL
  x.mds <- get_markers_from_ordered_sequence(x, lg = seq_along(x$maps), "mds")
  x.genome <- get_markers_from_ordered_sequence(x, lg = seq_along(x$maps), type = "genome")
  for(i in 1:length(x.mds)){
    a <- match(x.genome[[i]],x.mds[[i]])
    d <- rbind(d, data.frame(lg = names(x$maps)[i], x = seq_along(a), y = a))
  }
  # Modify this part to include your custom color palette
  num_colors <- length(unique(d$lg))  # Number of unique levels in lg
  color_palette <- mp_pal(num_colors) # Get the right number of colors from your palette

  p <- ggplot2::ggplot(d, ggplot2::aes(x, y)) +
    ggplot2::geom_point(alpha = alpha, ggplot2::aes(colour = lg), size = size) +
    ggplot2::scale_colour_manual(values = color_palette) + # Apply custom color palette
    ggplot2::facet_wrap(~lg, nrow = floor(sqrt(length(x$maps))), scales="free") +
    ggplot2::labs(subtitle = "Linkage group", x = "Genome position", y = "MDS position") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.position = "none", plot.subtitle = ggplot2::element_text(hjust = 0.5))

  # Display the plot
  print(p)
}

#' @export
plot.mappoly2.order.comparison <- function(x, ...){
  op <- par(xpd = TRUE)
  on.exit(par(op))
  pal <- c("#56B4E9","#E69F00")
  w <- unlist(x$maps, recursive = FALSE)
  max.dist <- max(unlist(x$maps))
  max.dist <- max.dist * 1.1
  plot(0,
       xlim = c(0, max.dist),
       ylim = c(0,length(w)+1),
       type = "n", axes = FALSE,
       xlab = "Map position (cM)",
       ylab = "Linkage groups")
  axis(1)
  for(i in 1:length(w)){
    plot_one_map(w[[i]] + max.dist * 0.1, i = i, horiz = TRUE, col = pal[1.5+((-1)^i)/2])
  }
  axis(2, at = 1:length(w), labels = names(w), lwd = 0, las = 2, hadj = 0)
  legend("topright",
         c("MDS",
           "Genome"),
         col = pal,
         pch = c(15, 15), horiz = TRUE, inset=c(0,-0.1))
  par(xpd = FALSE)
}

#' Plot Multiple Genetic Maps in a Grid Layout
#'
#' This function creates a visual representation of genetic maps from different biparental populations.
#' Each panel in the grid corresponds to a linkage group or chromosome. The maps are displayed
#' with the marker positions in centimorgans on the Y axis, and the distinct biparental populations
#' on the X axis. A color gradient from light to dark blue indicates the frequency of markers shared
#' across populations, providing insight into genetic similarities and differences.
#'
#' @param x A list of 'mappoly2.sequence' objects representing the genetic data of biparental populations.
#'
#' @return A ggplot object representing the genetic maps in a grid layout. Each panel shows a
#'         linkage group or chromosome with marker positions and shared frequency across populations.
#' @importFrom ggplot2 ggplot aes geom_point xlab ylab theme scale_color_manual theme margin guides guide_legend theme_dark
#' @importFrom dplyr count left_join vars
#' @export
plot_multi_map <- function(x){
  # Ensure all elements in 'x' are of class 'mappoly2.sequence'
  assert_that(all(sapply(x, function(x) is.mappoly2.sequence(x))),
              msg = "all elements in 'x' must be of class 'mappoly2.sequence'")

  # Construct names for each biparental population
  names(x) <- sapply(x, function(x) paste0(x$data$name.p1 , "x", x$data$name.p2))

  # Create a matrix to identify available maps for each population
  map.mat <- sapply(x, function(x)  sapply(x$maps, function(x) !is.null(x[["genome"]][["p1p2"]])), simplify = "array")
  map.mat <- matrix(map.mat, nrow = length(x[[1]]$maps), dimnames = list(names(x[[1]]$maps), names(x)))

  # Check for populations without maps and throw an error if any are found
  y <- apply(map.mat, 1, function(x) all(!x))
  if(any(y))
    stop("at least one population should have a map for group(s): ", paste(names(y)[y], collapse = " "))

  # Initialize a variable to store map data
  maps <- NULL

  # Loop through each population and extract map data
  for(i in 1:length(x)){
    w <- x[[i]]
    for(j in rownames(map.mat)[map.mat[,i]]){
      maps <- rbind(maps, data.frame(F1 = names(x)[i],
                                     LG = j,
                                     mrk.names = rownames(w$maps[[j]]$genome$p1p2$hmm.phase[[1]]$p1),
                                     pos = cumsum(c(0, imf_h(w$maps[[j]]$genome$p1p2$hmm.phase[[1]]$rf)))))
    }
  }

  mrk.names <- n <- pos <- F1 <- category <- LG <- NULL

  # Count the occurrences of each marker name in 'mrk.names'
  marker_counts <- maps %>%
    count(mrk.names) %>%
    mutate(category = as.character(n))

  # Join the counts back to the original data frame
  maps_with_counts <- maps %>%
    left_join(marker_counts, by = "mrk.names")

  # Generate a color set for each unique count category using viridis
  unique_counts <- sort(unique(marker_counts$category))
  colors <- rev(viridis::viridis(length(unique_counts)))
  names(colors) <- unique_counts

  # Update the plotting code
  lo <- optimal_layout(nrow(map.mat))
  ggplot(maps_with_counts, aes(x = pos, y = F1, group = as.factor(F1), color = category)) +
    geom_point(shape = 108, size = 5, show.legend = TRUE) +
    facet_wrap(vars(LG), nrow = lo[1], ncol = lo[2], scales = "free_x") +
    xlab("Position (cM)") +
    ylab("Biparental Maps") +
    scale_color_manual(values = colors, name = "Marker\nFrequency\nAcross\nPopulations") +  # Set the legend title
    theme(legend.title = element_text(face = "bold")) +  # Optionally make the legend title bold
    guides(color = guide_legend(override.aes = list(shape = 15))) + # Square shape for legend keys
    theme_dark()
}

#' Plot Consensus Map
#'
#' This function plots consensus genetic maps along with individual population maps.
#' It allows the option to plot only the consensus map or include individual population maps.
#'
#' @param x A list containing elements of 'mappoly2.prepared.integrated.data' class,
#'          which includes both individual maps and consensus map data.
#' @param only.consensus Logical, if TRUE, only the consensus map is plotted;
#'          if FALSE, individual population maps are included (default is FALSE).
#' @param col The color used for plotting the consensus map markers when
#' 'only.consensus' is TRUE (default is "lightgray").
#' @param ... Additional arguments, not used in this method
#'
#' @return A ggplot object representing the plotted genetic map(s).
#'
#' @importFrom ggplot2 ggplot geom_point facet_wrap xlab ylab scale_color_manual theme
#' @importFrom dplyr filter count mutate left_join
#' @export
plot.mappoly2.consensus.map <- function(x, only.consensus = FALSE, col = "lightgray", ...){
  z <- x$consensus.map
  x <- x$individual.maps
  # Ensure all elements in 'x' are of class 'mappoly2.sequence'
  assert_that(all(sapply(x, function(x) is.mappoly2.sequence(x))),
              msg = "all elements in 'x' must be of class 'mappoly2.sequence'")

  # Construct names for each biparental population
  names(x) <- sapply(x, function(x) paste0(x$data$name.p1 , "x", x$data$name.p2))

  # Create a matrix to identify available maps for each population
  map.mat <- sapply(x, function(x)  sapply(x$maps, function(x) !is.null(x[["genome"]][["p1p2"]])), simplify = "array")
  map.mat <- matrix(map.mat, nrow = length(x[[1]]$maps), dimnames = list(names(x[[1]]$maps), names(x)))

  # Check for populations without maps and throw an error if any are found
  y <- apply(map.mat, 1, function(x) all(!x))
  if(any(y))
    stop("at least one population should have a map for group(s): ", paste(names(y)[y], collapse = " "))

  if(only.consensus){
    if(length(col) == 1)
      col <- rep(col, nrow(map.mat))
    max.dist <- max(sapply(z, function(x) sum(imf_h(x$rf))))
    plot(0,
         xlim = c(0, max.dist),
         ylim = c(0,nrow(map.mat)+1),
         type = "n", axes = FALSE,
         xlab = "Map position (cM)",
         ylab = "",
         main = "Consensus Map")
    axis(1)
    axis(2, at = 1:nrow(map.mat), labels = rownames(map.mat), lwd = 0, las = 2)
    for(i in 1:nrow(map.mat)){
      d <- cumsum(c(0, imf_h(z[[rownames(map.mat)[i]]]$rf)))
      plot_one_map(d, i = i, horiz = TRUE, col = col[i])
    }
  }
  else{
    # Initialize a variable to store map data
    maps1 <- NULL

    # Loop through each population and extract map data
    for(i in 1:length(x)){
      w <- x[[i]]
      for(j in rownames(map.mat)[map.mat[,i]]){
        maps1 <- rbind(maps1, data.frame(POP = names(x)[i],
                                         LG = j,
                                         mrk.names = rownames(w$maps[[j]]$genome$p1p2$hmm.phase[[1]]$p1),
                                         pos = cumsum(c(0, imf_h(w$maps[[j]]$genome$p1p2$hmm.phase[[1]]$rf)))))
      }
    }
    # Append the consensus map data to the maps data frame
    maps2 <- NULL
    ## Consensus map
    for(j in names(z)){
      maps2 <- rbind(maps2, data.frame(POP = "Consensus",
                                       LG = j,
                                       mrk.names = rownames(z[[j]]$ph$PH[[1]]),
                                       pos = cumsum(c(0, imf_h(z[[j]]$rf)))))
    }
    maps <- rbind(maps1, maps2)

    # Adjust the factor levels of POP so that Consensus is first
    maps$POP <- factor(maps$POP, levels = c("Consensus", unique(maps$POP[maps$POP != "Consensus"])))

    # Continue with your existing code
    mrk.names <- n <- pos <- POP <- category <- LG <- NULL

    # Count the occurrences of each marker name in 'mrk.names', excluding 'Consensus'
    marker_counts <- maps %>%
      filter(POP != "Consensus") %>%
      count(mrk.names) %>%
      mutate(category = as.factor(n))

    # Join the counts back to the original data frame
    maps_with_counts <- maps %>%
      left_join(marker_counts, by = "mrk.names")

    # Generate a color set for each unique count category using viridis, excluding 'Consensus'
    unique_counts <- sort(unique(marker_counts$category))
    colors <- viridis::viridis(length(unique_counts), direction = -1)
    names(colors) <- unique_counts

    # Update the plotting code
    lo <- optimal_layout(nrow(map.mat))
    ggplot(maps_with_counts, aes(x = pos, y = POP, group = as.factor(POP), color = category)) +
      geom_point(shape = 108, size = 5, show.legend = TRUE) +
      facet_wrap(vars(LG), nrow = lo[1], ncol = lo[2]) +
      xlab("Position (cM)") +
      ylab("") +
      scale_color_manual(values = colors, name = "Marker\nFrequency\nAcross\nPopulations") +
      theme(legend.title = element_text(face = "bold")) +
      guides(color = guide_legend(override.aes = list(shape = 15))) +  # Square shape for legend keys
      theme_dark()

  }
}

#' Plot Haplotype Probabilities for Genetic Maps
#'
#' This function plots haplotype probabilities for specified linkage groups in a genetic mapping dataset.
#' It supports visualization of haplotype probabilities for individual or multiple linkage groups.
#'
#' @param x An object representing genetic mapping data, typically of a class storing genetic information.
#' @param lg Optional vector specifying the linkage group indices to plot.
#'           If NULL, the first linkage group is processed.
#' @param ind Index or name of the individual for which haplotypes are to be plotted.
#'            Should be a single numeric index or a character name.
#' @param type Character vector indicating the type of map to process, either "mds" or "genome".
#' @param parent Character vector specifying the parent or parents to be considered
#'               in the haplotype probability visualization. Options are "p1p2", "p1", and "p2".
#'
#' @return A ggplot object representing the plotted haplotype probabilities.
#'
#' @importFrom ggplot2 ggplot geom_density facet_grid theme_minimal ylab xlab scale_color_identity scale_fill_identity
#' @importFrom reshape2 melt
#' @importFrom assertthat assert_that
#' @export
plot_haplotypes <- function(x,
                            lg = NULL,
                            ind = 1,
                            type = c("mds", "genome"),
                            parent = c("p1p2", "p1", "p2")) {
  assert_that(is.mappoly2.sequence(x))
  if(is.null(lg))
    lg <- 1:length(x$maps)

  type <- match.arg(type)
  parent <- match.arg(parent)

  y <- parse_lg_and_type(x,lg,type)

  for(i in 1:length(y$lg))
    assert_that(is.haplotype.sequence(x, y$lg[i], y$type, parent),
                msg = "Requested haplotype probabilities were not computed")

  ### Asserting haplotype
  v <- detect_comp_haplotype(x)
  u <- apply(v[parent,,],1,all)
  h <- names(u)[1:2][!u[1:2]]
  if(length(h) == 1){
    assert_that(u[y$type], msg = paste(h, "order has not been computed for", parent))
  } else {
    assert_that(u[y$type], msg = paste(h[1], "and", h[2],"orders have not been computed for", parent))
  }

  if(is.character(ind) & length(ind) == 1){
    ind.num <- match(ind, x$data$screened.data$ind.names)
  } else if(is.numeric(ind) & length(ind) == 1){
    ind.num <- ind
    ind <- x$data$screened.data$ind.names[ind]
  } else {
    stop("indivdual is not valid")
  }
  if(is.na(ind.num))
    stop("individual not found")

  all_parents <- as.factor(c(x$data$name.p1, x$data$name.p2))

  hap <- lapply(x$maps[y$lg], function(map_item){
    M <- as.matrix(map_item[[y$type]][[parent]]$hmm.phase[[1]]$haploprob)
    ind.num <- match(ind, x$data$screened.data$ind.names)
    id <- which(ind.num==M[,2])
    hom_df <- as.data.frame(M[id,])
    map.pos <- cumsum(imf_h(c(0, map_item[[y$type]][[parent]]$hmm.phase[[1]]$rf)))
    map.pos <- rep(map.pos, each = length(id))
    # Adding column names for clarity
    colnames(hom_df) <- c("parent", "ind", "homolog", paste0("V", 4:ncol(hom_df)))

    # Using melt from reshape2 to transform the data
    # id.vars are the columns that identify each row (here, parent, ind, homolog)
    # measure.vars are the columns that contain the measurements to be reshaped (in this case, the rest of the columns)
    long_df <- reshape2::melt(hom_df, id.vars = c("parent", "homolog"), measure.vars = colnames(hom_df)[4:ncol(hom_df)])

    # Creating the final dataframe with required columns
    result_df <- data.frame(parent = as.factor(all_parents[long_df$parent]),
                            homolog = as.factor(paste0("h", long_df$homolog)),
                            prob = long_df$value,
                            map.pos = map.pos,
                            parent_pos = rep(c(rep("p1", x$data$ploidy.p1), rep("p2", x$data$ploidy.p2)),
                                             length(map.pos)/(x$data$ploidy.p1 + x$data$ploidy.p2)))
  })
  result_df <- reshape2::melt(hap,
                              id.vars = c("parent", "homolog",  "prob",  "map.pos", "parent_pos"),
                              value.name = "LG")

  # Create a new column for colors
  result_df$color <- NA

  # Assign colors to each homolog within each parent
  for (p in levels(result_df$parent)) {
    homolog_levels <- levels(result_df$homolog[result_df$parent == p])
    palette <- get_palette(all_parents, p, length(homolog_levels))
    for (h in homolog_levels) {
      result_df$color[result_df$parent == p & result_df$homolog == h] <- palette[which(homolog_levels == h)]
    }
  }
  result_df$parent <- paste0(result_df$parent, "-", result_df$parent_pos)
  color <- map.pos <- prob <- NULL
  if(length(y$lg) == 1){
    p <- ggplot(result_df, aes(x = map.pos, y = prob, color = color, fill = color)) +
      geom_density(stat = "identity", alpha = 0.7) +
      ggplot2::ggtitle(paste(ind, "   LG", y$lg)) +
      facet_grid(parent + homolog  ~ .) +
      theme_minimal() +
      ylab("Homologs Probability") +
      xlab("Map Position") +
      scale_color_identity() +
      scale_fill_identity()
  }else{
    p <- ggplot(result_df, aes(x = map.pos, y = prob, color = color, fill = color)) +
      ggplot2::geom_density(stat = "identity", alpha = 0.7, position = "stack") +
      ggplot2::ggtitle(ind) +
      #ggplot2::facet_grid(rows = ggplot2::vars(LG)) +
      facet_grid(L1 + parent  ~ .) +
      theme_minimal() +
      ylab("Homologs Probability") +
      xlab("Map Position") +
      scale_color_identity() +
      scale_fill_identity()
  }
  return(p)
}


#' Plot Consensus Map with Homolog Probabilities
#'
#' This function generates a density plot for homolog probabilities across map positions
#' for a specified individual within a given linkage group. It uses haplotype probability data
#' from a consensus map to create the plot, with different colors representing different parents.
#'
#' @param x A list containing consensus map data, specifically `consensus.map`, `ph`, `pedigree`, and `haploprob` structures.
#' @param lg An integer specifying the linkage group to plot. Defaults to 1.
#' @param ind An identifier for the individual to plot. This can be either a numeric index or a character string corresponding to the individual's name. If it is a character string, the function matches it with the rownames of the pedigree. Defaults to 1.
#' @param ... Additional arguments to pass on to the plotting function.
#'
#' @return A ggplot object representing the density plot of homolog probabilities for the specified individual across map positions.
#'
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_density ggtitle facet_grid theme_minimal ylab xlab scale_color_identity scale_fill_identity
#' @importFrom reshape2 melt
plot_consensus_haplo <- function(x,
                                 lg = 1,
                                 ind = 1,
                                 ...){
  ##FIXME allow for multiple lg's
  pedigree <- x$consensus.map[[lg]]$ph$pedigree
  all_parents <- as.factor(names(x$consensus.map[[lg]]$ph$PH))

  M <- as.matrix(x$consensus.map[[lg]]$haploprob)

  if(is.character(ind) & length(ind) == 1){
    ind.num <- match(ind, rownames(pedigree))
  } else if(is.numeric(ind) & length(ind) == 1){
    ind.num <- ind
    ind <- rownames(pedigree)[ind]
  } else {
    stop("indivdual is not valid")
  }
  if(is.na(ind.num))
    stop("individual not found")

  id <- which(ind.num==M[,2])
  hom_df <- as.data.frame(M[id,])
  map.pos <- cumsum(imf_h(c(0, x$consensus.map[[lg]]$rf)))
  map.pos <- rep(map.pos, each = length(id))

  # Adding column names for clarity
  colnames(hom_df) <- c("parent", "ind", "homolog", paste0("V", 4:ncol(hom_df)))

  # Using melt from reshape2 to transform the data
  # id.vars are the columns that identify each row (here, parent, ind, homolog)
  # measure.vars are the columns that contain the measurements to be reshaped (in this case, the rest of the columns)
  long_df <- reshape2::melt(hom_df, id.vars = c("parent", "homolog"), measure.vars = colnames(hom_df)[4:ncol(hom_df)])

  # Creating the final dataframe with required columns
  result_df <- data.frame(parent = as.factor(names(x$consensus.map[[lg]]$ph$PH))[long_df$parent],
                          homolog = as.factor(paste0("h", long_df$homolog)),
                          prob = long_df$value,
                          map.pos = map.pos,
                          parent_pos = rep(c(rep("p1", pedigree[ind.num,3]),
                                             rep("p2", pedigree[ind.num,4])),
                                           length(map.pos)/sum(pedigree[ind.num,3:4])))

  # Create a new column for colors
  result_df$color <- NA
  color <- prob <- NULL
  # Assign colors to each homolog within each parent
  for (p in levels(result_df$parent)) {
    homolog_levels <- levels(result_df$homolog[result_df$parent == p])
    palette <- get_palette(all_parents, parent = p, length(homolog_levels))
    for (h in homolog_levels) {
      result_df$color[result_df$parent == p & result_df$homolog == h] <- palette[which(homolog_levels == h)]
    }
  }
  result_df$parent <- paste0(result_df$parent, "-", result_df$parent_pos)
  # Plotting
  p <- ggplot(result_df, aes(x = map.pos, y = prob, color = color, fill = color)) +
    geom_density(stat = "identity", alpha = 0.7) +
    ggplot2::ggtitle(paste(ind, "   LG", lg)) +
    facet_grid(parent + homolog  ~ .) +
    theme_minimal() +
    ylab("Homologs Probability") +
    xlab("Map Position") +
    scale_color_identity() +
    scale_fill_identity()# Use the actual colors assigned in the color column

  return(p)
}


#' Plot Shared Markers in Genetic Maps
#'
#' This function visualizes shared markers across multiple genetic maps
#' in a 'mappoly2.prepared.integrated.data' object. It uses Euler diagrams to represent the intersection of markers.
#'
#' @param x An object of class 'mappoly2.prepared.integrated.data'.
#'          It should contain individual maps from which shared markers are to be identified.
#'
#' @return An Euler diagram plot showing shared markers across the individual maps in `x`.
#'
#' @details The function extracts markers from each map within the provided 'mappoly2.prepared.integrated.data'
#'          object and then plots an Euler diagram to show the intersections (shared markers) among these maps.
#'          This visualization helps in understanding the overlap of genetic information across different maps.
#'
#' @importFrom eulerr euler
#' @importFrom graphics plot
#' @export
plot_shared_markers <- function(x){
  assert_that(inherits(x, "mappoly2.prepared.integrated.data"))
  set_list <- sapply(x$individual.maps, function(z)
    unlist(sapply(z$maps, function(x)
      rownames(x$genome$p1p2$hmm.phase[[1]]$p1)),
      use.names = FALSE))
  colors <- mp_pal(length(set_list))
  plot(euler(set_list, shape = "circle"), quantities = TRUE, fills = colors)
}




