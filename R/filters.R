#' Filter out redundant markers
#'
#' Filter out markers with identical dosage information for all individuals.
#'
#' @param input.seq an object of class \code{mappoly.sequence}
#' @return An object of class \code{mappoly2.seq}
#'
#' @export filter_redundant
#' @import graphics
filter_redundant <- function(input.seq)
{
  x<-input.seq$data$geno.dose[input.seq$mrk.names, ]
  dim(x)
  id <- duplicated(x, dimnames = TRUE)
  dat.unique <- input.seq$data$geno.dose[!id, ]
  if(nrow(x) == nrow(dat.unique))
    return(input.seq)
  dat.duplicated <- input.seq$data$geno.dose[id, , drop = FALSE]
  n1 <- apply(dat.unique, 1, paste, collapse = "")
  n2 <- apply(dat.duplicated, 1, paste, collapse = "")
  structure(list(mrk.names = rownames(dat.unique),
                 redundant = data.frame(kept = rownames(dat.unique)[match(n2,n1)],
                                        removed = rownames(dat.duplicated)),
                 data = input.seq$data),
            class = "mappoly2.sequence")
}

#' Filter missing genotypes
#'
#' Excludes markers or individuals based on their proportion of missing data
#'
#' @param input.data an object of class \code{mappoly.data}
#'
#' @param type one of the following options:
#' \itemize{
#'   \item \code{'marker'}{filter out markers based on their percentage of missing data (default)}
#'   \item \code{'individual'}{filter out individuals based on their percentage of missing data}
#' }
#' Please notice that removing individuals with certain amount of data can change some marker parameters
#' (such as depth), and can also change the estimated genotypes for other individuals.
#' So be careful when removing individuals.
#'
#' @param filter.thres maximum percentage of missing data (default = 0.2)
#'
#' @param inter if \code{TRUE}, expects user-input to proceed with filtering
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
#' @importFrom graphics axis
filter_missing <- function(input.data,
                           type = c("marker", "individual"),
                           filter.thres = 0.2,
                           inter = TRUE) {
  assert_that(is.mappoly2.data(input.data))
  type <- match.arg(type)
  switch(type,
         marker = filter_missing_mrk(input.data,
                                     filter.thres = filter.thres,
                                     inter = inter),
         individual = filter_missing_ind(input.data,
                                         filter.thres = filter.thres,
                                         inter = inter)
  )
}

#' Filter markers based on missing genotypes
#'
#' @param input.data an object of class \code{"mappoly.data"}
#' @param filter.thres maximum percentage of missing data
#' @param inter if \code{TRUE}, expects user-input to proceed with filtering
#' @keywords internal
filter_missing_mrk <- function(input.data, filter.thres = 0.2, inter = TRUE) {
  op <- par(pty = "s")
  on.exit(par(op))

  process_filter <- function(filter.thres) {
    na.num <- apply(input.data$geno.dose, 1, function(x) sum(is.na(x)))
    perc.na <- na.num / input.data$n.ind
    mrks.id <- which(perc.na <= filter.thres)

    if (length(mrks.id) == input.data$n.mrk) {
      return(input.data)
    }

    out.dat <- subset_data(input.data, select.mrk = names(mrks.id))
    return(out.dat)
  }

  if (interactive() && inter) {
    ANSWER <- "flag"

    while (!grepl("^(y|yes)?$", ANSWER, ignore.case = TRUE)) {
      na.num <- apply(input.data$geno.dose, 1, function(x) sum(is.na(x)))
      perc.na <- na.num / input.data$n.ind
      x <- sort(perc.na)

      plot(x,
           xlab = "markers",
           ylab = "frequency of missing data",
           col = ifelse(x <= filter.thres, 4, 2),
           pch = ifelse(x <= filter.thres, 1, 4))
      abline(h = filter.thres, lty = 2)

      legend("topleft",
             c(paste0("Filtered out: ", sum(perc.na > filter.thres)),
               paste0("Included: ", sum(perc.na <= filter.thres))),
             col = c(2, 4),
             pch = c(4, 1))

      ANSWER <- readline("Enter 'Y/n' to proceed or update the filter threshold: ")

      if (grepl("^(n|no)$", ANSWER, ignore.case = TRUE)) {
        stop("Stop function.")
      } else if (!grepl("^(y|yes)?$", ANSWER, ignore.case = TRUE)) {
        filter.thres <- as.numeric(ANSWER)
      }
    }
  }
  par(pty="m")
  return(process_filter(filter.thres))
}


#' Filter individuals based on missing genotypes
#'
#' @param input.data an object of class \code{"mappoly.data"}
#' @param filter.thres maximum percentage of missing data
#' @param inter if \code{TRUE}, expects user-input to proceed with filtering
#' @keywords internal
#' @importFrom graphics axis
filter_missing_ind <- function(input.data, filter.thres = 0.2, inter = TRUE) {
  op <- par(pty = "s")
  on.exit(par(op))

  process_filter <- function(filter.thres) {
    na.num <- apply(input.data$geno.dose, 2, function(x) sum(is.na(x)))
    perc.na <- na.num / input.data$n.mrk
    ind.id <- which(perc.na <= filter.thres)

    if (length(ind.id) == input.data$n.ind) {
      return(input.data)
    }

    out.dat <- subset_data(input.data, select.ind = names(ind.id))
    return(out.dat)
  }

  if (interactive() && inter) {
    ANSWER <- "flag"

    while (!grepl("^(y|yes)?$", ANSWER, ignore.case = TRUE)) {
      na.num <- apply(input.data$geno.dose, 2, function(x) sum(is.na(x)))
      perc.na <- na.num / input.data$n.mrk
      x <- sort(perc.na)

      plot(x,
           xlab = "offspring",
           ylab = "frequency of missing data",
           col = ifelse(x <= filter.thres, 4, 2),
           pch = ifelse(x <= filter.thres, 1, 4))
      abline(h = filter.thres, lty = 2)

      legend("topleft",
             c(paste0("Filtered out: ", sum(perc.na > filter.thres)),
               paste0("Included: ", sum(perc.na <= filter.thres))),
             col = c(2, 4),
             pch = c(4, 1))

      ANSWER <- readline("Enter 'Y/n' to proceed or update the filter threshold: ")

      if (grepl("^(n|no)$", ANSWER, ignore.case = TRUE)) {
        stop("You decided to stop the function.")
      } else if (!grepl("^(y|yes)?$", ANSWER, ignore.case = TRUE)) {
        filter.thres <- as.numeric(ANSWER)
      }
    }
  }
  par(pty="m")
  return(process_filter(filter.thres))
}


#' Filter Markers Based on Chi-Square Test
#'
#' This function filters markers based on the p-values obtained from a chi-square test.
#' The chi-square test assumes that the markers follow the expected segregation patterns under
#' Mendelian inheritance, random chromosome bivalent pairing, and no double reduction.
#'
#' @param input.obj The name of the input object of class \code{mappoly.data}.
#'
#' @param chisq.pval.thres The p-value threshold used for the chi-square tests. By default,
#' this is the Bonferroni approximation with a global alpha of 0.05, i.e., 0.05/n.mrk.
#'
#' @param inter If set to TRUE (default), this function will plot the distorted and
#' non-distorted markers.
#'
#' @return An object of class \code{mappoly.chitest.seq} that includes the following components:
#' \item{keep}{markers that follow Mendelian segregation pattern}
#' \item{exclude}{markers with distorted segregation}
#' \item{chisq.pval.thres}{threshold p-value used for chi-square tests}
#' \item{data.name}{name of the input dataset used to perform the chi-square tests}
#'
#' @author Marcelo Mollinari
#' \email{mmollin@ncsu.edu}
#'
#' @importFrom graphics axis
#' @export
filter_segregation <- function(input.obj, chisq.pval.thres = NULL, inter = TRUE){
  # Set plot options
  op <- par(pty="s")
  on.exit(par(op))

  # Extract chi-square p-values and number of markers
  if(is.mappoly2.data(input.obj)){
    chisq.val <- input.obj$chisq.pval
    n.mrk <- input.obj$n.mrk
  } else if (is.mappoly2.sequence(input.obj)){
    chisq.val <- input.obj$data$chisq.pval[input.obj$mrk.names]
    n.mrk <- length(input.obj$mrk.names)
  } else {
    stop(deparse(substitute(input.obj)),
         " is not an object of class 'mappoly2.data' or 'mappoly2.sequence'")
  }

  # Set threshold for chi-square p-values using Bonferroni approximation if not specified
  if(is.null(chisq.pval.thres))
    chisq.pval.thres <- 0.05/n.mrk

  # Prompt user for input if interactive and inter = TRUE
  ANSWER <- "flag"
  if(interactive() && inter) {
    while(substr(ANSWER, 1, 1) != "y" && substr(ANSWER, 1, 1) != "yes" && substr(ANSWER, 1, 1) != "Y" && ANSWER  != "") {
      x <- log10(sort(chisq.val, decreasing = TRUE))
      th <- log10(chisq.pval.thres)
      plot(x,
           xlab = "markers",
           ylab = bquote(log[10](P)),
           col = ifelse(x <= th, 2, 4),
           pch =ifelse(x <= th, 4, 1))
      abline(h = th, lty = 2)
      f <- paste0("Filtered out: ", sum(x < th))
      i <- paste0("Included: ", sum(x >= th))
      legend("bottomleft",  c(f, i) , col = c(2,4), pch = c(4,1))
      ANSWER <- readline("Enter 'Y/n' to proceed or update the p value threshold: ")
      if(substr(ANSWER, 1, 1)  ==  "n" | substr(ANSWER, 1, 1)  ==  "no" | substr(ANSWER, 1, 1)  ==  "N") {
        stop("You decided to stop the function.")
      }
      if(substr(ANSWER, 1, 1) != "y" && substr(ANSWER, 1, 1) != "yes" && substr(ANSWER, 1, 1) != "Y" && ANSWER  != "") {
        chisq.pval.thres  <- as.numeric(ANSWER)
      }
    }
  }

  # Identify markers that meet threshold for chi-square p-value and return filtered object
  keep <- names(which(chisq.val >= chisq.pval.thres))
  par(pty="m")
  return(make_sequence(input.obj, keep))
}


#' Filter out individuals
#'
#' This function removes individuals from the input dataset, either by specifying
#' them manually or by using interactive kinship analysis.
#'
#' @param input.data The name of the input object (class \code{mappoly.data}).
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
filter_individuals <- function(input.data,
                               ind.to.remove = NULL,
                               inter = TRUE,
                               verbose = TRUE){
  assert_that(is.mappoly2.data(input.data))
  if(input.data$ploidy.p1 != input.data$ploidy.p2)
    stop("'filter_individuals' cannot be executed\n  on progenies with odd ploidy levels.")
  op <- par(pty="s")
  on.exit(par(op))
  D <- t(input.data$geno.dose)
  D <- rbind(input.data$dosage.p1, input.data$dosage.p2, D)
  rownames(D)[1:2] <- c("P1", "P2")
  G  <- AGHmatrix::Gmatrix(D, method = "VanRaden",ploidy = input.data$ploidy.p1/2 + input.data$ploidy.p2/2)
  x <- G[1,]
  y <- G[2,]
  df <- data.frame(x = x, y = y, type = c(2, 2, rep(4, length(x)-2)))
  plot(df[,1:2], col = df$type, pch = 19,
       xlab = paste0("relationships between the offspring and ",input.data$name.p1),
       ylab = paste0("relationships between the offspring and ",input.data$name.p2))
  abline(c(0,1), lty = 2)
  abline(c(-0.4,1), lty = 2, col = "gray")
  abline(c(0.4,1), lty = 2, col = "gray")
  legend("topright",  c("Parents", "Offspring") , col = c(2,4), pch = 19)
  if(!is.null(ind.to.remove)){
    out.dat <- subset_data(input.data,  select.ind = setdiff(input.data$ind.names, ind.to.remove))
    return(out.dat)
  }
  if(interactive() && inter)
  {
 #   if (!require(gatepoints))
 #     stop("Please install package 'gatepoints' to proceed")
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
        return(input.data)
      }
      out.dat <- subset_data(input.data,  select.ind = ind.to.include)
      return(out.dat)
    } else{
      warning("No individuals removed. Returning original data set.")
      par(pty="m")
      return(input.data)
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
