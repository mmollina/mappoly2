#' @export
#' @importFrom graphics barplot layout mtext image legend
#' @importFrom grDevices colorRampPalette
plot.mappoly2 <- function(x, type = c("original", "filtered"), ...)
{
  type <- match.arg(type)
  if(type == "original"){
    assert_that(inherits(x, "data"))
      plot_data(x, type = "original")
  } else if(type == "filtered"){
    assert_that(inherits(x, "screened"))
      plot_data(x, type = "filtered", mrk.id = x$screened.data$mrk.names,
                ind.id = x$screened.data$ind.names)

  }
}

plot_data <- function(x, type, mrk.id = NULL, ind.id = NULL, ...){
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
  type2 <- apply(d.temp, 1, function(x,ploidy.p1, ploidy.p2) paste0(sort(abs(abs(as.numeric(x)-(ploidy.p1/2))-(ploidy.p2/2))), collapse = ""),
                ploidy.p1 = x$data$ploidy.p1, ploidy.p2 = x$data$ploidy.p2)
  type2.names <- names(table(type2))
  mrk.dist <- as.numeric(freq)
  names(mrk.dist) <- apply(d.temp, 1 , paste, collapse = "-")
  layout(matrix(c(1,1,1,2,3,3,6,4,5), 3, 3), widths = c(1.2,3,.5), heights = c(1.5,2,3))
  barplot(mrk.dist, las = 2, #col = pal[match(type, type.names)],
          xlab = "Number of markers",
          ylab = "Dosage combination", horiz = TRUE)
  pval <- x$data$screening$markers[,"chisq.pval"]
  names(pval) <- rownames(x$data$screening$markers)
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
  if(!is.null(x$data$redundant) & type != "filtered"){
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

