
#' Filter Out Redundant Markers
#'
#' This function filters out redundant markers from a `mappoly2.data` object. It identifies and removes duplicated genetic markers based on their genotypic information, thus reducing the data to unique markers. It also updates the `redundant` component of the object with information about which markers were kept and which were removed.
#'
#' @param x An object of class \code{mappoly2.data}. It is expected to have a component named `geno.dose`, which is a matrix containing the dosage for each marker (rows) for each individual (columns).
#'
#' @return The modified `mappoly2.data` object with redundant markers removed. The return value includes:
#'   - The original data with redundant markers removed.
#'   - The `redundant` component of the object updated to reflect changes. This component is a data frame listing the markers that were kept (`kept`) and the corresponding markers that were removed due to redundancy (`removed`).
#'
#' @details This function first checks if the input is a valid `mappoly2.data` object using \code{assert_that(is.mappoly2.data(x))}. It then identifies redundant markers by checking for duplicates in the `geno.dose` matrix. If no redundant markers are found, the function returns the original object. If redundant markers are identified, they are removed, and the object's `redundant` data frame is updated accordingly.
#' @keywords internal
filter_redundant <- function(x)
{
  assert_that(is.mappoly2.data(x))
  id <- duplicated(x$geno.dose, dimnames = TRUE)
  dat.unique <- x$geno.dose[!id, ]
  if(nrow(x$geno.dose) == nrow(dat.unique)){
    if(!is.data.frame(x$redundant))
      x$redundant <- 0
    else{
      w <- x$redundant
      x$redundant <- w[w$kept%in%rownames(x$geno.dose),]
    }
    return(x)
  }
  dat.duplicated <- x$geno.dose[id, , drop = FALSE]
  n1 <- apply(dat.unique, 1, paste, collapse = "")
  n2 <- apply(dat.duplicated, 1, paste, collapse = "")
  w <- data.frame(kept = rownames(dat.unique)[match(n2,n1)],
                  removed = rownames(dat.duplicated))
  x <- subset_data(x, select.mrk = setdiff(x$mrk.names,
                                           w$removed))
  w <- unique(rbind(x$redundant, w))
  x$redundant <- w[w$kept%in%rownames(x$geno.dose),]
  return(x)
}

#' Filter Genetic Data Based on Quality Metrics
#'
#' This function filters genetic data in a `mappoly2.data` object based on various quality control metrics. It applies thresholds for missing data rates, chi-squared p-values, and read depth to markers and individuals, and optionally plots the screening process.
#'
#' @param x A `mappoly2.data` object containing genetic data.
#' @param mrk.thresh A numeric threshold for the missing data rate in markers (default is 0.10).
#' @param ind.thresh A numeric threshold for the missing data rate in individuals (default is 0.10).
#' @param chisq.pval.thresh A numeric threshold for chi-squared test p-values in markers (default is NULL, which sets the threshold using a Bonferroni approximation).
#' @param read.depth.thresh A numeric vector with two values indicating the lower and upper bounds for acceptable read depths in markers (default is c(5, 1000)).
#' @param plot.screening Logical, if TRUE (default), plots are generated to visually represent the screening process.
#'
#' @return Returns the input `mappoly2.data` object with additional components:
#'   - `screened.data`: A list containing the thresholds used for selection and the names of markers and individuals that met the specified criteria.
#'   - Class attribute `screened` is also appended to the object.
#'
#' @details The function first validates the input object, then applies the specified thresholds to filter out markers and individuals based on missing data rates, chi-squared p-values, and read depths. It updates the object with the results of this filtering and optionally generates plots to visualize the data before and after filtering.
#'
#' @examples
#'   filtered_data <- filter_data(B2721)
#'
#' @importFrom graphics axis
#' @export
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
  id <- .get_mrk_ind_from_QAQC(x$QAQC.values,
                               miss.mrk.thresh = mrk.thresh,
                               miss.ind.thresh = ind.thresh,
                               chisq.pval.thresh = chisq.pval.thresh,
                               read.depth.thresh = read.depth.thresh)
  x$screened.data <- id
  x$screened.data$thresholds <- c(x$screened.data$thresholds,
                                  LOD.ph = 0,
                                  LOD.rf = 0,
                                  rf = 0.5,
                                  prob.lower = 0,
                                  prob.upper = 1)
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
#' @importFrom AGHmatrix Gmatrix
#' @importFrom gatepoints fhs
#' @export
filter_individuals <- function(x,
                               ind.to.remove = NULL,
                               inter = TRUE,
                               type = c("Gmat", "PCA"),
                               verbose = TRUE){
  assert_that(is.mappoly2.data(x))
  if(x$ploidy.p1 != x$ploidy.p2)
    stop("'filter_individuals' cannot be executed\n  on progenies with odd ploidy levels.")
  type <- match.arg(type)
  op <- par(pty="s")
  on.exit(par(op))
  D <- t(x$geno.dose)
  if(has.mappoly2.screened(x)){
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
      id <- .get_mrk_ind_from_QAQC(x$QAQC.values,
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
        id <- .get_mrk_ind_from_QAQC(x$QAQC.values,
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

#' Remove Markers Not Meeting LOD and Recombination Fraction Criteria
#'
#' This function removes markers from a `mappoly2.data` or `mappoly2.sequence` object
#' that do not meet specified LOD (logarithm of odds) and recombination fraction criteria.
#' It is designed to filter out markers that are unlikely to be linked or show strong
#' evidence of linkage across an entire linkage group, which might indicate false positives.
#'
#' @param x An object of class \code{mappoly2.data} or \code{mappoly2.sequence}.
#'
#' @param thresh.LOD.ph LOD score threshold for linkage phase configuration.
#' Typically set to eliminate markers that are unlinked to the analyzed set (default = 5).
#'
#' @param thresh.LOD.rf LOD score threshold for recombination fraction (default = 5).
#'
#' @param thresh.rf Recombination fraction threshold (default = 0.15).
#'
#' @param probs A numeric vector indicating the probability corresponding to the filtering quantiles (default = c(0.05, 1)).
#'
#' @param lg A vector of linkage groups to be processed.
#' If NULL (default), all groups are considered.
#'
#' @param type A character vector specifying the method for linkage group analysis.
#' Options are "mds" for multidimensional scaling and "genome" for genomic analysis.
#' This parameter only has an effect if `lg` is not NULL.
#'
#' @param diag.markers A vector specifying a window of marker pairs to consider.
#' If NULL (default), all markers are considered.
#'
#' @param mrk.order Marker order vector.
#' This parameter is only used if `diag.markers` is not NULL.
#'
#' @param diagnostic.plot Logical; if \code{TRUE}, generates a diagnostic plot (default is TRUE).
#'
#' @param breaks The number of cells for the histogram in the diagnostic plot (default = 100).
#'
#' @return A filtered object of the same class as the input (`mappoly2.data` or `mappoly2.sequence`).
#'
#' @details The function first checks the type of the input object and applies the
#' relevant filtering criteria based on LOD scores and recombination fractions.
#' It optionally produces a diagnostic plot to visualize the filtering process.
#'
#' @examples
#' \dontrun{
#'   # Assuming `my_data` is a valid mappoly2.data or mappoly2.sequence object
#'   filtered_data <- rf_filter(my_data, thresh.LOD.ph = 5, thresh.LOD.rf = 5)
#' }
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu} with updates by Gabriel Gesteira, \email{gdesiqu@ncsu.edu}
#' @export
#' @importFrom ggplot2 ggplot geom_histogram aes scale_fill_manual xlab ggtitle
#' @importFrom graphics hist par text points
rf_filter <- function(x,
                      thresh.LOD.ph = 5,
                      thresh.LOD.rf = 5,
                      thresh.rf = 0.15,
                      probs = c(0.05, 1),
                      lg = NULL,
                      type = c("mds", "genome"),
                      diag.markers = NULL,
                      mrk.order = NULL,
                      diagnostic.plot = TRUE,
                      breaks = 100){
  if(is.mappoly2.data(x)){
    assert_that(has.mappoly2.rf(x))
    return(init_rf_filter(x,
                          thresh.LOD.ph,
                          thresh.LOD.rf,
                          thresh.rf,
                          probs,
                          diag.markers,
                          mrk.order,
                          diagnostic.plot,
                          breaks = breaks))
  } else if(is.mappoly2.sequence(x)){
    lg.temp <- c(1 : length(x$maps))
    if(is.null(lg)){
      lg <-lg.temp
      diagnostic.plot <- FALSE
    }
    assert_that(all(lg %in% lg.temp), msg = "Provide a valid group set")
    for(i in lg){
      x <- rf_filter_per_group(x,
                               i,
                               type,
                               thresh.LOD.ph,
                               thresh.LOD.rf,
                               thresh.rf,
                               probs,
                               diag.markers,
                               diagnostic.plot,
                               breaks)
    }
    return(x)
  }
}

#' Initialize Recombination Fraction Filtering
#'
#' Internal function to initialize filtering based on recombination fractions, LOD scores for linkage phase and recombination fraction.
#'
#' @param x An object of class \code{pairwise.rf}.
#' @param thresh.LOD.ph LOD threshold for linkage phase (default = 5).
#' @param thresh.LOD.rf LOD threshold for recombination fraction (default = 5).
#' @param thresh.rf Recombination fraction threshold (default = 0.15).
#' @param probs Probability range for filtering (default = c(0.05, 1)).
#' @param diag.markers Window of marker pairs to consider (default = NULL).
#' @param mrk.order Order of markers (used if diag.markers is not NULL).
#' @param diagnostic.plot Boolean to control the generation of a diagnostic plot (default = TRUE).
#' @param breaks Number of breaks for histogram in diagnostic plot (default = 100).
#' @return A modified \code{pairwise.rf} object with filtered data.
#' @keywords internal
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
  assert_that(inherits(x, "pairwise.rf"))
  probs <- range(probs)
  ## Getting filtered rf matrix
  M <-filter_rf_matrix(x,
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
  x$screened.data$thresholds$LOD.ph = thresh.LOD.ph
  x$screened.data$thresholds$LOD.rf = thresh.LOD.rf
  x$screened.data$thresholds$rf = thresh.rf
  x$screened.data$thresholds$prob.lower = probs[1]
  x$screened.data$thresholds$prob.upper = probs[2]
  x$screened.data$mrk.names <- intersect(x$screened.data$mrk.names, ids)
  return(x)
}

#' Filter Recombination Fractions for Specific Linkage Group
#'
#' Internal function to filter recombination fractions for a specific linkage group in a \code{mappoly2.sequence} object.
#'
#' @param x A \code{mappoly2.sequence} object.
#' @param lg Linkage group to filter.
#' @param type Method for linkage group analysis ("mds" or "genome").
#' @param thresh.LOD.ph LOD threshold for linkage phase (default = 5).
#' @param thresh.LOD.rf LOD threshold for recombination fraction (default = 5).
#' @param thresh.rf Recombination fraction threshold (default = 0.15).
#' @param probs Probability range for filtering (default = c(0.05, 1)).
#' @param diag.markers Window of marker pairs to consider (default = NULL).
#' @param diagnostic.plot Boolean to control the generation of a diagnostic plot (default = TRUE).
#' @param breaks Number of breaks for histogram in diagnostic plot (default = 100).
#' @return A modified \code{mappoly2.sequence} object with filtered linkage group.
#' @keywords internal
rf_filter_per_group <- function(x,
                                lg,
                                type = c("mds", "genome"),
                                thresh.LOD.ph = 5,
                                thresh.LOD.rf = 5,
                                thresh.rf = 0.15,
                                probs = c(0.05, 1),
                                diag.markers = NULL,
                                diagnostic.plot = TRUE,
                                breaks = 100)
{
  y <- parse_lg_and_type(x,lg,type)
  assert_that(length(y$lg) ==1 & is.numeric(lg))
  probs <- range(probs)
  ## Getting filtered rf matrix
  mrk.names <- x$maps[[lg]][[y$type]]$mkr.names
  M <-filter_rf_matrix(x$data,
                       type = "rf",
                       thresh.LOD.ph = thresh.LOD.ph,
                       thresh.LOD.rf = thresh.LOD.rf,
                       thresh.rf = thresh.rf,
                       mrk.names = mrk.names)
  if(!is.null(diag.markers))
    M[abs(col(M) - row(M)) > diag.markers] <- NA
  z <- apply(M, 1, function(x) sum(!is.na(x)))
  w <- hist(z, breaks = breaks, plot = FALSE)
  th <- quantile(z, probs = probs)
  rem <- c(which(z < th[1]), which(z > th[2]))
  ids <- names(which(z >= th[1] & z <= th[2]))
  value <- marker <- NULL
  if(diagnostic.plot){
    d <- rbind(data.frame(marker = "original", value = z),
               data.frame(marker = "filtered", value = z[ids]))
    p <- ggplot2::ggplot(d, ggplot2::aes(value)) +
      ggplot2::geom_histogram(ggplot2::aes(fill = marker),
                              alpha = 0.4, position = "identity", binwidth = diff(w$mids)[1]) +
      ggplot2::scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
      ggplot2::ggtitle( paste0("Filtering probs: [", probs[1], " : ", probs[2], "] - Non NA values by row in rf matrix - b width: ", diff(w$mids)[1])) +
      ggplot2::xlab(paste0("Non 'NA' values at LOD.ph = ", thresh.LOD.ph, ", LOD.rf = ", thresh.LOD.rf, ", and thresh.rf = ", thresh.rf))
    print(p)
  }
  x$maps[[lg]][[y$type]]$mkr.names <- intersect(x$maps[[lg]][[y$type]]$mkr.names, ids)
  return(x)
}

#' Filter Recombination Fraction Matrix
#'
#' Internal function to filter a recombination fraction matrix based on LOD thresholds and recombination fraction limits.
#'
#' @param x A data object containing pairwise recombination fractions.
#' @param type Type of matrix to filter ("rf" for recombination fraction, "sh" for SH values).
#' @param thresh.LOD.ph LOD threshold for linkage phase (default = 0).
#' @param thresh.LOD.rf LOD threshold for recombination fraction (default = 0).
#' @param thresh.rf Recombination fraction threshold (default = 0.5).
#' @param mrk.names Marker names to include in the filtering (default = NULL, all markers).
#' @return A filtered recombination fraction matrix or a list of filtered SH matrices.
#' @keywords internal
filter_rf_matrix <- function(x,
                             type = c("rf", "sh"),
                             thresh.LOD.ph = 0,
                             thresh.LOD.rf = 0,
                             thresh.rf = 0.5,
                             mrk.names = NULL){
  type <- match.arg(type)
  if(is.null(mrk.names))
    mrk.names <- colnames(x$pairwise.rf$rec.mat)
  lod.ph.mat <- x$pairwise.rf$lod.ph.mat[mrk.names,mrk.names]
  lod.mat <- x$pairwise.rf$lod.mat[mrk.names,mrk.names]
  rec.mat <- x$pairwise.rf$rec.mat[mrk.names,mrk.names]
  id1 <- abs(lod.ph.mat) < thresh.LOD.ph
  id2 <- abs(lod.mat) < thresh.LOD.rf
  id3 <- rec.mat > thresh.rf
  if(type == "rf"){
    if(thresh.LOD.ph > 0) rec.mat[id1] <- NA
    if(thresh.LOD.rf > 0) rec.mat[id2] <- NA
    if(thresh.rf < 0.5) rec.mat[id3] <- NA
    return(rec.mat)
  } else if(type == "sh"){
    sh.p1 <- x$pairwise.rf$Sh.p1[mrk.names,mrk.names]
    sh.p2 <- x$pairwise.rf$Sh.p2[mrk.names,mrk.names]
    if(thresh.LOD.ph > 0) sh.p2[id1] <- sh.p1[id1] <- NA
    if(thresh.LOD.rf > 0) sh.p2[id2] <- sh.p1[id2] <- NA
    if(thresh.rf < 0.5)  sh.p2[id3] <- sh.p1[id3] <- NA
    return(list(Sh.p1 = sh.p1,
                Sh.p2 = sh.p2))
  }
}

#' Edit sequence ordered by reference genome positions
#' comparing to another set order
#'
#' @param input.seq object of class mappoly2.sequence with alternative order (not genomic order)
#' @param group linkage group id
#' @param invert vector of marker names to be inverted
#' @param remove vector of marker names to be removed
#'
#'  
#' @export
edit_order <- function(input.seq, group = 1,invert = NULL, remove = NULL){
  
  if (!inherits(input.seq, "mappoly2.sequence")) {
    stop(deparse(substitute(input.seq)), " is not an object of class 'mappoly2.sequence'")
  }
  
  if (is.null(input.seq$maps$lg1$mds$order) | is.null(input.seq$maps$lg1$genome$order)) {
    stop("Run `order_sequence` with type = `mds` and type = `genome` before editing your sequence")
  }
  
  d <- lg <- y <- NULL
  x.mds <- get_markers_from_ordered_sequence(input.seq, lg = group, 
                                             "mds")
  x.genome <- get_markers_from_ordered_sequence(input.seq, lg = group, 
                                                type = "genome")
  for (i in 1:length(x.mds)) {
    a <- match(x.genome[[i]], x.mds[[i]])
    d <- rbind(d, data.frame(x = seq_along(a), 
                             y = a))
  }
  
  rownames(d) <- x.genome[[1]]
  plot(d$x, d$y, xlab="Genome position", ylab = "MDS position")
  
  inverted <- removed <- vector()
  if(!is.null(invert) | !is.null(remove)){
    if(!is.null(invert)){
      inverted <- c(inverted, as.vector(invert))
      repl <- d[rev(match(as.vector(invert),rownames(d))),]
      d[match(as.vector(invert),rownames(d)),2] <- repl[,2]
      rownames(d)[match(as.vector(invert), rownames(d))] <- rownames(repl) 
    }
    if(!is.null(remove)){
      removed <- c(removed, as.vector(remove))
      d <- d[-match(remove, rownames(d)),]
    }
    plot(d$x, d$y, xlab="Genome position", ylab = "MDS position")
  } else {
    cat("Mark at least three points on the plot and press `Esc` to continue.")
    if(interactive()){
      ANSWER <- "Y"
      while(substr(ANSWER, 1, 1)  ==  "y" | substr(ANSWER, 1, 1)  ==  "yes" | substr(ANSWER, 1, 1)  ==  "Y" | ANSWER  == ""){
        plot(d$x, d$y, xlab="Genome position", ylab = "MDS position")
        mks.to.remove <- gatepoints::fhs(d, mark = TRUE)
        if(length(which(rownames(d) %in% mks.to.remove)) > 0){
          ANSWER2 <- readline("Enter 'invert/remove' to proceed with the edition: ")
          if(ANSWER2 == "invert"){
            inverted <- c(inverted, as.vector(mks.to.remove))
            repl <- d[rev(match(as.vector(mks.to.remove),rownames(d))),]
            d[match(as.vector(mks.to.remove),rownames(d)),2] <- repl[,2]
            rownames(d)[match(as.vector(mks.to.remove), rownames(d))] <- rownames(repl) 
          } else {
            removed <- c(removed, as.vector(mks.to.remove))
            d <- d[-match(mks.to.remove, rownames(d)),]
          }
        }
        ANSWER <- readline("Enter 'Y/n' to proceed with interactive edition or quit: ")
      }
      plot(d$x, d$y, xlab="Genome position", ylab = "MDS position")
    }
  }
  
  custom_ord <- input.seq$maps[[group]]$genome$order[match(rownames(d), rownames(input.seq$maps[[group]]$genome$order)),]
  input.seq$maps[[group]]$custom$order <- custom_ord
  return(input.seq)
}