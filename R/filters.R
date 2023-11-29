.setQAQC <- function(id.mrk, id.ind,
                     miss.mrk = NA,
                     miss.ind = rep(NA, length(id.ind)),
                     chisq.pval){
  list(markers = data.frame(miss = miss.mrk,
                            chisq.pval = chisq.pval,
                            read.depth = NA,
                            row.names = id.mrk),
       individuals = data.frame(miss = miss.ind,
                                full.sib = NA,
                                row.names = id.ind))
}

.get_mrk_ind_from_QAQC <- function(x,
                                   miss.mrk.thresh = +Inf,
                                   miss.ind.thresh = +Inf,
                                   chisq.pval.thresh = -Inf,
                                   read.depth.thresh = c(0, +Inf)){

  if(any(is.na(x$markers[,"read.depth"]))){
    id.mrk <- x$markers[,"miss"] < miss.mrk.thresh &
      x$markers[,"chisq.pval"] > chisq.pval.thresh
  } else{
    id.mrk <- x$markers[,"miss"] < miss.mrk.thresh &
      x$markers[,"chisq.pval"] > chisq.pval.thresh &
      x$markers[,"read.depth"] > read.depth.thresh[1] &
      x$markers[,"read.depth"] < read.depth.thresh[2]
  }
  if(any(is.na(x$individuals[,"full.sib"])))
    id.ind <- x$individuals[,"miss"] < miss.ind.thresh
  else
    id.ind <- x$individuals[,"miss"] < miss.ind.thresh &
      x$individuals[,"full.sib"]
  return(list(thresholds = list(miss.mrk = miss.mrk.thresh,
                                miss.ind = miss.ind.thresh,
                                chisq.pval = chisq.pval.thresh,
                                read.depth = read.depth.thresh),
              mrk.names = rownames(x$markers)[id.mrk],
              ind.names = rownames(x$individuals)[id.ind]))
}

#' Filter out markers with redundant information
filter_redundant <- function(x)
{
  id <- duplicated(x$geno.dose, dimnames = TRUE)
  dat.unique <- x$geno.dose[!id, ]
  if(nrow(x$geno.dose) == nrow(dat.unique))
    return(NA)
  dat.duplicated <- x$geno.dose[id, , drop = FALSE]
  n1 <- apply(dat.unique, 1, paste, collapse = "")
  n2 <- apply(dat.duplicated, 1, paste, collapse = "")
  return(data.frame(kept = rownames(dat.unique)[match(n2,n1)],
                    removed = rownames(dat.duplicated)))
}

#' @export
#' @importFrom graphics axis
filter_data <- function(x,
                        mrk.thresh = 0.10,
                        ind.thresh = 0.10,
                        chisq.pval.thresh = NULL,
                        read.depth.thresh = c(5,1000),
                        plot.screening = TRUE) {
  assert_that(inherits(x, "mappoly2.data"))
  op <- par(pty = "s", mfrow = c(2,2), mar = c(4,3,3,2))
  on.exit(par(op))
  chisq.val <- x$QAQC.values$markers$chisq.pval
  # Set threshold for chi-square p-values using Bonferroni approximation if not specified
  if(is.null(chisq.pval.thresh))
    chisq.pval.thresh <- 0.05/length(chisq.val)
  id <- mappoly2:::.get_mrk_ind_from_QAQC(x$QAQC.values,
                                          miss.mrk.thresh = mrk.thresh,
                                          miss.ind.thresh = ind.thresh,
                                          chisq.pval.thresh = chisq.pval.thresh,
                                          read.depth.thresh = read.depth.thresh)
  x$screened.data <- id
  class(x) <- c(class(x), "screened")
  pal <- c("#56B4E9","#E69F00")
  if (plot.screening) {
    ####Missing markers ####
    z <- sort(x$QAQC.values$markers$miss)
    rg <- range(z)
    if(rg[2] < .1) rg[2] <- .1
    plot(z,
         xlab = "markers",
         ylab = "",
         col = ifelse(z <= mrk.thresh, pal[1], pal[2]),
         pch = ifelse(z <= mrk.thresh, 1, 4),
         main = "Markers", xlim = c(ceiling(-length(z)*.05) , ceiling(length(z)*1.05)),
         ylim = rg)
    mtext("frequency of missing data", 2, line = 2, cex= .75)
    abline(h = mrk.thresh, lty = 2)
    legend("topleft",
           c(paste0("Filtered out: ", sum(z > mrk.thresh)),
             paste0("Included: ", sum(z <= mrk.thresh))),
           col = rev(pal),
           pch = c(4, 1))
    ####Missing individuals ####
    z <- sort(x$QAQC.values$individuals$miss)
    rg <- range(z)
    if(rg[2] < .1) rg[2] <- .1
    plot(z,
         xlab = "individuals",
         ylab = "",
         col = ifelse(z <= ind.thresh, pal[1], pal[2]),
         pch = ifelse(z <= ind.thresh, 1, 4),
         main = "Individuals",
         xlim = c(ceiling(-length(z)*.05) , ceiling(length(z)*1.05)),
         ylim = rg)
    mtext("frequency of missing data", 2, line = 2, cex= .75)
    abline(h = ind.thresh, lty = 2)
    legend("topleft",
           c(paste0("Filtered out: ", sum(z > ind.thresh)),
             paste0("Included: ", sum(z <= ind.thresh))),
           col = rev(pal),
           pch = c(4, 1))
    #### Chi-square test ####
    w <- log10(sort(x$QAQC.values$markers$chisq.pval, decreasing = TRUE))
    th <- log10(chisq.pval.thresh)
    plot(w,
         xlab = "markers",
         ylab = "",
         col = ifelse(w <= th, pal[2], pal[1]),
         pch =ifelse(w <= th, 4, 1), main = "Segregation", xlim = c(ceiling(-length(w)*.05) , ceiling(length(w)*1.05)))
    abline(h = th, lty = 2)
    mtext(bquote(log[10](P)), 2, line = 2, cex= .75)
    f <- paste0("Filtered out: ", sum(w < th))
    i <- paste0("Included: ", sum(w >= th))
    legend("bottomleft",  c(f, i) , col = rev(pal), pch = c(4,1))
    #### Read Depth ####
    if(all(!is.na(x$QAQC.values$markers$read.depth))){
      hist_info <- hist(x$QAQC.values$markers$read.depth,
                        main = "Read depth", xlab = "number of reads", col = pal[1])
      lower_tail <- read.depth.thresh[1]
      upper_tail <- read.depth.thresh[2]
      if(hist_info$breaks[2] <  lower_tail)
      {
        # Add colored rectangles for the tails
        rect(hist_info$breaks[hist_info$breaks < lower_tail],
             0,
             hist_info$breaks[which(hist_info$breaks < lower_tail) + 1],
             hist_info$counts[hist_info$breaks < lower_tail],
             col=pal[2])
      }
      if(rev(hist_info$breaks)[2] >  upper_tail)
      {
        rect(hist_info$breaks[hist_info$breaks >= upper_tail],
             0,
             hist_info$breaks[which(hist_info$breaks >= upper_tail) + 1],
             hist_info$counts[hist_info$breaks >= upper_tail],
             col=pal[2])
      }
      abline(v = read.depth.thresh, lty = 2)
    } else {
      plot(0, type = "n", axes = FALSE, xlab = "", ylab = "")
    }
    par(pty="m")
  }
  return(x)
}

#' Filter out individuals
#'
#' This function removes individuals from the input dataset, either by specifying
#' them manually or by using interactive kinship analysis.
#'
#' @param x The name of the input object (class \code{mappoly.data}).
#'
#' @param ind.to.remove A character vector containing the names of the individuals
#' to be removed. If \code{NULL}, the function opens an interactive graphic to
#' allow for individual selection.
#'
#' @param inter If \code{TRUE}, the function expects user input to proceed with
#' filtering.
#'
#' @param type A character string specifying the procedure to be used for
#'             detecting outlier offspring. Options include "Gmat", which
#'             utilizes the genomic kinship matrix, and "PCA", which employs
#'             principal component analysis on the dosage matrix.
#'
#' @param verbose If \code{TRUE} (default), the function shows the list of filtered
#'  out individuals.
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @importFrom stats prcomp
#' @export
filter_individuals <- function(x,
                               ind.to.remove = NULL,
                               inter = TRUE,
                               type = c("Gmat", "PCA"),
                               verbose = TRUE){
  assert_that(mappoly2:::is.mappoly2.data(x))
  if(x$ploidy.p1 != x$ploidy.p2)
    stop("'filter_individuals' cannot be executed\n  on progenies with odd ploidy levels.")
  type <- match.arg(type)
  op <- par(pty="s")
  on.exit(par(op))
  D <- t(x$geno.dose)
  if(is.mappol2.screened(x)){
    D <- D[x$screened.data$ind.names, x$screened.data$mrk.names]
    D <- rbind(x$dosage.p1[x$screened.data$mrk.names],
               x$dosage.p2[x$screened.data$mrk.names],
               D)
  } else {
    D <- rbind(x$dosage.p1,
               x$dosage.p2,
               D)
  }
  rownames(D)[1:2] <- c(x$name.p1, x$name.p2)
  if(type == "Gmat"){
    G  <- AGHmatrix::Gmatrix(D, method = "VanRaden",ploidy = x$ploidy.p1/2 + x$ploidy.p2/2)
    y1 <- G[1,]
    y2 <- G[2,]
    df <- data.frame(x = y1, y = y2, type = c(2, 2, rep(4, length(y1)-2)))
    plot(df[,1:2], col = df$type, pch = 19,
         xlab = paste0("relationships between the offspring and ",x$name.p1),
         ylab = paste0("relationships between the offspring and ",x$name.p2))
    abline(c(0,1), lty = 2)
    abline(c(-0.4,1), lty = 2, col = "gray")
    abline(c(0.4,1), lty = 2, col = "gray")
    legend("topright",  c("Parents", "Offspring") , col = c(2,4), pch = 19)
  }
  else{
    row_means <- rowMeans(D, na.rm = TRUE)
    for (i in 1:nrow(D))
      D[i, is.na(D[i, ])] <- row_means[i]
    pc <- prcomp(D)
    x <- pc$x[,"PC1"]
    y <- pc$x[,"PC2"]
    a <- diff(range(x))*0.05
    b <- diff(range(y))*0.05
    df <- data.frame(x = x, y = y, type = c(2, 2, rep(4, length(x)-2)))
    plot(df[,1:2], col = df$type, pch = 19,
         xlab = "PC1",
         ylab = "PC2",
         xlim = c(min(x)-a, max(x)+a),
         ylim = c(min(y)-b, max(y)+b))
    points(df[1:2,1:2], col = 2, pch = 19)
    legend("bottomleft",  c("Parents", "Offspring") , col = c(2,4), pch = 19)
  }
  if(!is.null(ind.to.remove)){
    full.sib <- !x$ind.names%in%ind.to.remove
    x$QAQC.values$individuals[,"full.sib"] <- !rownames(x$QAQC.values$individuals)%in%ind.to.remove
    if(inherits(x, "screened")){
      id <- mappoly2:::.get_mrk_ind_from_QAQC(x$QAQC.values,
                                              miss.mrk.thresh = x$screened.data$thresholds$miss.mrk,
                                              miss.ind.thresh = x$screened.data$thresholds$miss.ind,
                                              chisq.pval.thresh = x$screened.data$thresholds$chisq.pval,
                                              read.depth.thresh = x$screened.data$thresholds$read.depth)
    }
    x$screened.data <- id
    return(x)
  }
  if(interactive() && inter)
  {
    ANSWER <- readline("Enter 'Y/n' to proceed with interactive filtering or quit: ")
    if(substr(ANSWER, 1, 1)  ==  "y" | substr(ANSWER, 1, 1)  ==  "yes" | substr(ANSWER, 1, 1)  ==  "Y" | ANSWER  == "")
    {
      ind.to.remove <- gatepoints::fhs(df, mark = TRUE)
      ind.to.remove <- setdiff(ind.to.remove, c("P1", "P2"))
      ind.to.include <- setdiff(rownames(df)[-c(1:2)], ind.to.remove)
      if(verbose){
        cat("Removing individual(s): \n")
        print(ind.to.remove)
        cat("...\n")
      }
      if(length(ind.to.remove) == 0){
        warning("No individuals removed. Returning original data set.")
        par(pty="m")
        return(x)
      }
      full.sib <- !x$ind.names%in%ind.to.remove
      x$QAQC.values$individuals[,"full.sib"] <- !rownames(x$QAQC.values$individuals)%in%ind.to.remove
      if(inherits(x, "screened")){
        id <- mappoly2:::.get_mrk_ind_from_QAQC(x$QAQC.values,
                                                miss.mrk.thresh = x$screened.data$thresholds$miss.mrk,
                                                miss.ind.thresh = x$screened.data$thresholds$miss.ind,
                                                chisq.pval.thresh = x$screened.data$thresholds$chisq.pval,
                                                read.depth.thresh = x$screened.data$thresholds$read.depth)
      }
      x$screened.data <- id
      par(pty="m")
      return(x)
    } else{
      warning("No individuals removed. Returning original data set.")
      par(pty="m")
      return(x)
    }
  }
  par(pty="m")
}
#'  Remove markers that do not meet a LOD criteria
#'
#'  Remove markers that do not meet a LOD and recombination fraction
#'  criteria for at least a percentage of the pairwise marker
#'  combinations. It also removes markers with strong evidence of
#'  linkage across the whole linkage group (false positive).
#'
#' \code{thresh.LOD.ph} should be set in order to only select
#'     recombination fractions that have LOD scores associated to the
#'     linkage phase configuration higher than \code{thresh_LOD_ph}
#'     when compared to the second most likely linkage phase configuration.
#'     That action usually eliminates markers that are unlinked to the
#'     set of analyzed markers.
#'
#' @param input.twopt an object of class \code{mappoly.twopt}
#'
#' @param thresh.LOD.ph LOD score threshold for linkage phase configuration
#' (default = 5)
#'
#' @param thresh.LOD.rf LOD score threshold for recombination fraction
#' (default = 5)
#'
#' @param thresh.rf threshold for recombination fractions (default = 0.15)
#'
#' @param probs indicates the probability corresponding to the filtering
#' quantiles. (default = c(0.05, 1))
#'
#' @param diag.markers A window where marker pairs should be considered.
#'    If NULL (default), all markers are considered.
#'
#' @param mrk.order marker order. Only has effect if 'diag.markers' is not NULL
#'
#' @param ncpus number of parallel processes (i.e. cores) to spawn
#' (default = 1)
#'
#' @param diagnostic.plot if \code{TRUE} produces a diagnostic plot
#'
#' @param breaks number of cells for the histogram
#'
#' @return A filtered object of class \code{mappoly.sequence}.
#'
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu} with updates by Gabriel Gesteira, \email{gdesiqu@ncsu.edu}
#'
#' @export
#' @importFrom ggplot2 ggplot geom_histogram aes scale_fill_manual xlab ggtitle
#' @importFrom graphics hist
init_rf_filter <- function(x,
                           thresh.LOD.ph = 5,
                           thresh.LOD.rf = 5,
                           thresh.rf = 0.15,
                           probs = c(0.05, 1),
                           diag.markers = NULL,
                           mrk.order = NULL,
                           diagnostic.plot = TRUE,
                           breaks = 100)
{
  assert_that(inherits(x, "pairwise"))
  probs <- range(probs)
  ## Getting filtered rf matrix
  M <-mappoly2:::filter_rf_matrix(x,
                                  type = "rf",
                                  thresh.LOD.ph = thresh.LOD.ph,
                                  thresh.LOD.rf = thresh.LOD.rf,
                                  thresh.rf = thresh.rf)
  if(!is.null(mrk.order))
    M <- M[mrk.order, mrk.order]
  if(!is.null(diag.markers))
    M[abs(col(M) - row(M)) > diag.markers] <- NA
  z <- apply(M, 1, function(x) sum(!is.na(x)))
  w <- hist(z, breaks = breaks, plot = FALSE)
  th <- quantile(z, probs = probs)
  rem <- c(which(z < th[1]), which(z > th[2]))
  ids <- names(which(z >= th[1] & z <= th[2]))
  value <- type <- NULL
  if(diagnostic.plot){
    d <- rbind(data.frame(type = "original", value = z),
               data.frame(type = "filtered", value = z[ids]))
    p <- ggplot2::ggplot(d, ggplot2::aes(value)) +
      ggplot2::geom_histogram(ggplot2::aes(fill = type),
                              alpha = 0.4, position = "identity", binwidth = diff(w$mids)[1]) +
      ggplot2::scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
      ggplot2::ggtitle( paste0("Filtering probs: [", probs[1], " : ", probs[2], "] - Non NA values by row in rf matrix - b width: ", diff(w$mids)[1])) +
      ggplot2::xlab(paste0("Non 'NA' values at LOD.ph = ", thresh.LOD.ph, ", LOD.rf = ", thresh.LOD.rf, ", and thresh.rf = ", thresh.rf))
    print(p)
  }
  x$initial.screened.rf <- list(thresholds = c(thresh.LOD.ph = thresh.LOD.ph,
                                               thresh.LOD.rf = thresh.LOD.rf,
                                               thresh.rf = thresh.rf,
                                               prob.lower = probs[1],
                                               prob.upper = probs[2]),
                                mrk.names = ids)
  class(x) <- unique(c(class(x), "init_rf_screened"))
  return(x)
}

#' @export
rf_filter_per_group <- function(x,
                                gr,
                                thresh.LOD.ph = 5,
                                thresh.LOD.rf = 5,
                                thresh.rf = 0.15,
                                probs = c(0.05, 1),
                                diag.markers = NULL,
                                mrk.order = NULL,
                                diagnostic.plot = TRUE,
                                breaks = 100)
{
  assert_that(inherits(x, "ws"))
  mrk.names <- x$working.sequences[[gr]]$mrk.names
  probs <- range(probs)
  ## Getting filtered rf matrix
  M <- mappoly2:::filter_rf_matrix(x,
                                  type = "rf",
                                  thresh.LOD.ph = thresh.LOD.ph,
                                  thresh.LOD.rf = thresh.LOD.rf,
                                  thresh.rf = thresh.rf,
                                  mrk.names)

  if(!is.null(mrk.order)){
    assert_that(all(mrk.order%in%mrk.names))
    M <- M[mrk.order, mrk.order]
  }
  if(!is.null(diag.markers))
    M[abs(col(M) - row(M)) > diag.markers] <- NA
  z <- apply(M, 1, function(x) sum(!is.na(x)))
  w <- hist(z, breaks = breaks, plot = FALSE)
  th <- quantile(z, probs = probs)
  rem <- c(which(z < th[1]), which(z > th[2]))
  ids <- names(which(z >= th[1] & z <= th[2]))
  value <- type <- NULL
  if(diagnostic.plot){
    d <- rbind(data.frame(type = "original", value = z),
               data.frame(type = "filtered", value = z[ids]))
    p <- ggplot2::ggplot(d, ggplot2::aes(value)) +
      ggplot2::geom_histogram(ggplot2::aes(fill = type),
                              alpha = 0.4, position = "identity", binwidth = diff(w$mids)[1]) +
      ggplot2::scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
      ggplot2::ggtitle( paste0("Filtering probs: [", probs[1], " : ", probs[2], "] - Non NA values by row in rf matrix - b width: ", diff(w$mids)[1])) +
      ggplot2::xlab(paste0("Non 'NA' values at LOD.ph = ", thresh.LOD.ph, ", LOD.rf = ", thresh.LOD.rf, ", and thresh.rf = ", thresh.rf))
    print(p)
  }
  x$working.sequences[[gr]]$screened.rf <- list(thresholds = c(thresh.LOD.ph = thresh.LOD.ph,
                                                               thresh.LOD.rf = thresh.LOD.rf,
                                                               thresh.rf = thresh.rf,
                                                               prob.lower = probs[1],
                                                               prob.upper = probs[2]),
                                                mrk.names = ids)
  return(x)
}


filter_rf_matrix <- function(x,
                             type = c("rf", "sh"),
                             thresh.LOD.ph = 0,
                             thresh.LOD.rf = 0,
                             thresh.rf = 0.5,
                             mrk.names = NULL){
  type <- match.arg(type)
  if(is.null(mrk.names))
    mrk.names <- colnames(x$pairwise$rec.mat)
  lod.ph.mat <- x$pairwise$lod.ph.mat[mrk.names,mrk.names]
  lod.mat <- x$pairwise$lod.mat[mrk.names,mrk.names]
  rec.mat <- x$pairwise$rec.mat[mrk.names,mrk.names]
  id1 <- abs(lod.ph.mat) < thresh.LOD.ph
  id2 <- abs(lod.mat) < thresh.LOD.rf
  id3 <- rec.mat > thresh.rf
  if(type == "rf"){
    if(thresh.LOD.ph > 0) rec.mat[id1] <- NA
    if(thresh.LOD.rf > 0) rec.mat[id2] <- NA
    if(thresh.rf < 0.5) rec.mat[id3] <- NA
    return(rec.mat)
  } else if(type == "sh"){
    sh.p1 <- x$pairwise$Sh.p1[mrk.names,mrk.names]
    sh.p2 <- x$pairwise$Sh.p2[mrk.names,mrk.names]
    if(thresh.LOD.ph > 0) sh.p1[id1] <- sh.p1[id1] <- NA
    if(thresh.LOD.rf > 0) sh.p1[id2] <- sh.p1[id2] <- NA
    if(thresh.rf < 0.5)  sh.p1[id3] <- sh.p1[id3] <- NA
    return(list(Sh.p1 = sh.p1,
                Sh.p2 = sh.p2))
  }
}

