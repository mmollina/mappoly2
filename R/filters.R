.setScreeningClass <- function(id.mrk, id.ind,
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

.get_mrk_ind_from_Screening <- function(x,
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
  id <- duplicated(x$data$geno.dose, dimnames = TRUE)
  dat.unique <- x$data$geno.dose[!id, ]
  if(nrow(x$data$geno.dose) == nrow(dat.unique))
    return(NA)
  dat.duplicated <- x$data$geno.dose[id, , drop = FALSE]
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
  assert_that(inherits(x, "data"))
  assert_that(all(class(x)%in%c("mappoly2", "data", "screened")))
  op <- par(pty = "s", mfrow = c(2,2), mar = c(3,2,2,1))
  on.exit(par(op))
  chisq.val <- x$data$screening$markers$chisq.pval
  # Set threshold for chi-square p-values using Bonferroni approximation if not specified
  if(is.null(chisq.pval.thresh))
    chisq.pval.thresh <- 0.05/length(chisq.val)
  id <- mappoly2:::.get_mrk_ind_from_Screening(x$data$screening,
                                               miss.mrk.thresh = mrk.thresh,
                                               miss.ind.thresh = ind.thresh,
                                               chisq.pval.thresh = chisq.pval.thresh,
                                               read.depth.thresh = read.depth.thresh)
  x$screened.data <- id
  class(x) <- c(class(x), "screened")
  pal <- c("#56B4E9","#E69F00")
  if (plot.screening) {
    ####Missing markers ####
    z <- sort(x$data$screening$markers$miss)
    plot(z,
         xlab = "markers",
         ylab = "frequency of missing data",
         col = ifelse(z <= mrk.thresh, pal[1], pal[2]),
         pch = ifelse(z <= mrk.thresh, 1, 4),
         main = "Markers", xlim = c(ceiling(-length(z)*.05) , ceiling(length(z)*1.05)))
    abline(h = mrk.thresh, lty = 2)
    legend("topleft",
           c(paste0("Filtered out: ", sum(z > mrk.thresh)),
             paste0("Included: ", sum(z <= mrk.thresh))),
           col = rev(pal),
           pch = c(4, 1))
    ####Missing individuals ####
    z <- sort(x$data$screening$individuals$miss)
    plot(z,
         xlab = "individuals",
         ylab = "frequency of missing data",
         col = ifelse(z <= ind.thresh, pal[1], pal[2]),
         pch = ifelse(z <= ind.thresh, 1, 4),
         main = "Individuals", xlim = c(ceiling(-length(z)*.05) , ceiling(length(z)*1.05)))
    abline(h = ind.thresh, lty = 2)
    legend("topleft",
           c(paste0("Filtered out: ", sum(z > ind.thresh)),
             paste0("Included: ", sum(z <= ind.thresh))),
           col = rev(pal),
           pch = c(4, 1))
    #### Chi-square test ####
    w <- log10(sort(x$data$screening$markers$chisq.pval, decreasing = TRUE))
    th <- log10(chisq.pval.thresh)
    plot(w,
         xlab = "markers",
         ylab = bquote(log[10](P)),
         col = ifelse(w <= th, pal[2], pal[1]),
         pch =ifelse(w <= th, 4, 1), main = "Segregation", xlim = c(ceiling(-length(w)*.05) , ceiling(length(w)*1.05)))
    abline(h = th, lty = 2)
    f <- paste0("Filtered out: ", sum(w < th))
    i <- paste0("Included: ", sum(w >= th))
    legend("bottomleft",  c(f, i) , col = rev(pal), pch = c(4,1))
    #### Read Depth ####
    if(all(!is.na(x$data$screening$markers$read.depth))){
      hist_info <- hist(x$data$screening$markers$read.depth,
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
#' @param verbose If \code{TRUE} (default), the function shows the list of filtered
#'  out individuals.
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @export
filter_individuals <- function(x,
                               ind.to.remove = NULL,
                               inter = TRUE,
                               verbose = TRUE){
  assert_that(inherits(x, "data"))
  assert_that(all(class(x)%in%c("mappoly2", "data", "screened")))
  if(x$data$ploidy.p1 != x$data$ploidy.p2)
    stop("'filter_individuals' cannot be executed\n  on progenies with odd ploidy levels.")
  op <- par(pty="s")
  on.exit(par(op))
  D <- t(x$data$geno.dose)
  if(inherits(x, "screened")){
   D <- D[x$screened.data$ind.names, x$screened.data$mrk.names]
   D <- rbind(x$data$dosage.p1[x$screened.data$mrk.names],
              x$data$dosage.p2[x$screened.data$mrk.names],
              D)
  } else {
    D <- rbind(x$data$dosage.p1,
               x$data$dosage.p2,
               D)
  }
  rownames(D)[1:2] <- c(x$data$name.p1, x$data$name.p2)
  G  <- AGHmatrix::Gmatrix(D, method = "VanRaden",ploidy = x$data$ploidy.p1/2 + x$data$ploidy.p2/2)
  y1 <- G[1,]
  y2 <- G[2,]
  df <- data.frame(x = y1, y = y2, type = c(2, 2, rep(4, length(y1)-2)))
  plot(df[,1:2], col = df$type, pch = 19,
       xlab = paste0("relationships between the offspring and ",x$data$name.p1),
       ylab = paste0("relationships between the offspring and ",x$data$name.p2))
  abline(c(0,1), lty = 2)
  abline(c(-0.4,1), lty = 2, col = "gray")
  abline(c(0.4,1), lty = 2, col = "gray")
  legend("topright",  c("Parents", "Offspring") , col = c(2,4), pch = 19)
  if(!is.null(ind.to.remove)){
    full.sib <- !x$data$ind.names%in%ind.to.remove
    x$data$screening$individuals[,"full.sib"] <- !rownames(x$data$screening$individuals)%in%ind.to.remove
    if(inherits(x, "screened")){
    id <- mappoly2:::.get_mrk_ind_from_Screening(x$data$screening,
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
        return(x)
      }
      full.sib <- !x$data$ind.names%in%ind.to.remove
      x$data$screening$individuals[,"full.sib"] <- !rownames(x$data$screening$individuals)%in%ind.to.remove
      if(inherits(x, "screened")){
        id <- mappoly2:::.get_mrk_ind_from_Screening(x$data$screening,
                                                     miss.mrk.thresh = x$screened.data$thresholds$miss.mrk,
                                                     miss.ind.thresh = x$screened.data$thresholds$miss.ind,
                                                     chisq.pval.thresh = x$screened.data$thresholds$chisq.pval,
                                                     read.depth.thresh = x$screened.data$thresholds$read.depth)
      }
      x$screened.data <- id
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
#' @export rf_snp_filter
#' @importFrom ggplot2 ggplot geom_histogram aes scale_fill_manual xlab ggtitle
#' @importFrom graphics hist
rf_snp_filter <- function(input.twopt,
                          thresh.LOD.ph = 5,
                          thresh.LOD.rf = 5,
                          thresh.rf = 0.15,
                          probs = c(0.05, 1),
                          diag.markers = NULL,
                          mrk.order = NULL,
                          ncpus = 1L,
                          diagnostic.plot = TRUE,
                          breaks = 100)
{
  assert_that(is.mappoly2.twopt(input.twopt))
  probs <- range(probs)
  ## Getting filtered rf matrix
  rf_mat <-  rf_list_to_matrix(input.twopt = input.twopt, thresh.LOD.ph = thresh.LOD.ph,
                               thresh.LOD.rf = thresh.LOD.rf, thresh.rf = thresh.rf,
                               ncpus = ncpus, verbose = FALSE)
  M <- rf_mat$rec.mat
  if(!is.null(mrk.order))
    M <- M[mrk.order, mrk.order]
  if(!is.null(diag.markers))
    M[abs(col(M) - row(M)) > diag.markers] <- NA
  x <- apply(M, 1, function(x) sum(!is.na(x)))
  w <- hist(x, breaks = breaks, plot = FALSE)
  th <- quantile(x, probs = probs)
  rem <- c(which(x < th[1]), which(x > th[2]))
  ids <- names(which(x >= th[1] & x <= th[2]))
  value <- type <- NULL
  if(diagnostic.plot){
    d <- rbind(data.frame(type = "original", value = x),
               data.frame(type = "filtered", value = x[ids]))
    p <- ggplot2::ggplot(d, ggplot2::aes(value)) +
      ggplot2::geom_histogram(ggplot2::aes(fill = type),
                              alpha = 0.4, position = "identity", binwidth = diff(w$mids)[1]) +
      ggplot2::scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
      ggplot2::ggtitle( paste0("Filtering probs: [", probs[1], " : ", probs[2], "] - Non NA values by row in rf matrix - b width: ", diff(w$mids)[1])) +
      ggplot2::xlab(paste0("Non 'NA' values at LOD.ph = ", thresh.LOD.ph, ", LOD.rf = ", thresh.LOD.rf, ", and thresh.rf = ", thresh.rf))
    print(p)
  }
  ## Returning sequence object
  ch_filt <- make_sequence(input.obj = input.twopt$input.seq$data,
                           arg = ids)
  return(ch_filt)
}
