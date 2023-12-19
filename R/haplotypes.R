#' Calculate Haplotype probabilities using Hidden Markov Models
#'
#' This function calculates haplotypes for each linkage group in a genetic mapping dataset. It supports both sequential and parallel processing.
#'
#' @param x An object representing genetic mapping data.
#' @param lg Optional; a vector of linkage group indices to process. If NULL, all linkage groups in `x` are processed.
#' @param type The type of map to process, either "mds" or "genome".
#' @param phase.conf A configuration parameter for phase calculation.
#' @param verbose Logical; if TRUE, progress messages will be printed.
#' @param ncpus The number of CPU cores to use for parallel processing.
#' @return The input object `x` with haplotypes calculated for the specified linkage groups.
#' @importFrom parallel detectCores makeCluster parLapply stopCluster
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
calc_haplotypes <- function(x,
                            lg = NULL,
                            type = c("mds", "genome"),
                            parent = c("p1p2","p1","p2"),
                            phase.conf = "all",
                            verbose = TRUE,
                            ncpus = 1)
{
  y <- mappoly2:::parse_lg_and_type(x,lg,type)
  parent <- match.arg(parent)
  g <- x$data$geno.dose
  g[is.na(g)] <- -1
  ploidy.p1 <- x$data$ploidy.p1
  ploidy.p2 <- x$data$ploidy.p2
  dosage.p1 <- x$data$dosage.p1
  dosage.p2 <- x$data$dosage.p2
  ind.names <- x$data$screened.data$ind.names

  # Assessing multi-point map availability
  has.hmm.map <- sapply(x$maps[y$lg], function(x) sapply(x[[y$type]][3:5], function(x) !is.null(x$hmm.phase[[1]]$loglike)))

  # Checking for minimal phase information
  if(all(!has.hmm.map))
    stop("Provide a hmm estimated map.")

  ## Selecting markers based on input
  u <- t(has.hmm.map[parent, , drop = FALSE])
  v <- apply(u, 1, any)
  if(all(!v))
    stop(paste("Provide a pre-phased sequence for groups", paste(names(v), collapse = " ")))
  if(any(!v))
    warning(paste("Provide a pre-phased sequence for groups", paste(names(v[!v]), collapse = " ")))
  p <- apply(u[v,,drop = FALSE], 1, function(x) ifelse(x[1], "hmm.phase", NA))
  haplotypeData <- vector("list", length(p))
  names(haplotypeData) <- names(p)
  for(i in names(p)){
    ph <- x$maps[[i]][[y$type]][[parent]][[p[i]]]
    gtemp <- g[rownames(ph[[1]]$p1), ind.names]
    haplotypeData[[i]] <- list(g = gtemp,
                               ph = ph,
                               ploidy.p1 = ploidy.p1,
                               ploidy.p2 = ploidy.p2,
                               dosage.p1 = dosage.p1,
                               dosage.p2 = dosage.p2,
                               phase.conf = phase.conf,
                               parent.info = parent,
                               verbose = verbose)
  }
  haplotypeFunc <- function(data) {
    return(mappoly2:::calc_haplotypes_one(data$g, data$ph, data$ploidy.p1,
                                          data$ploidy.p2, dat$dosage.p1,
                                          dat$dosage.p2, data$phase.conf,
                                          parent.info = data$parent.info,
                                          verbose = data$verbose))
  }
  if(ncpus > 1) {
    cl <- makeCluster(ncpus)
    results <- parLapply(cl, haplotypeData, haplotypeFunc)
    stopCluster(cl)
  } else {
    results <- lapply(haplotypeData, haplotypeFunc)
  }
  for(i in names(p))
    x$maps[[i]][[y$type]][[parent]][["hmm.phase"]] <- results[[i]]
  return(x)
}

calc_haplotypes_one <- function(g,
                                ph,
                                ploidy.p1,
                                ploidy.p2,
                                dosage.p1,
                                dosage.p2,
                                phase.conf = "all",
                                parent.info = c("p1p2", "p1", "p2"),
                                verbose = TRUE){
  if(all(phase.conf == "all"))
    phase.conf <- 1:length(ph)
  assert_that(all(phase.conf%in%1:length(ph)),
              msg = "invalid phases specified in 'phase.conf'")
  n.ind <- ncol(g)
  mrk.id <- rownames(ph[[1]]$p1)
  g <- g[mrk.id,]
  if (parent.info == "p1p2"){ ###FIXME: include detect_parent_info
    for(i in phase.conf){
      cat("   Conf.", i,":")
      pedigree <- matrix(rep(c(1,2,ploidy.p1,ploidy.p2, 1),n.ind),
                         nrow = n.ind,
                         byrow = TRUE)
      w <- calc_haploprob_biallelic(PH = list(ph[[i]]$p1,
                                              ph[[i]]$p2),
                                    G = g,
                                    pedigree = pedigree,
                                    rf = ph[[i]]$rf,
                                    err = ph[[i]]$error)
      ph[[i]]$haploprob <- cbind(rep(c(rep(1, ploidy.p1), rep(2, ploidy.p2)), n.ind), w)
    }
    cat("Done with haplotype probabilities\n")
    return(ph)
  }
  else if(parent.info == "p1"){
    id <- which(ploidy.p2 == dosage.p2[mrk.id])
    g[id, ] <- g[id, ] - ploidy.p2/2
    for(i in phase.conf){
      cat("   Conf.", i,":")
      w <- calc_haploprob_biallelic_single(PH = ph[[i]]$p1,
                                           G = g,
                                           rf = ph[[i]]$rf,
                                           err = ph[[i]]$error)
      ph[[i]]$haploprob <- cbind(rep(1, ploidy.p1*n.ind), w)
    }
    cat("Done with haplotype probabilities\n")
    return(ph)
  }
  else if(parent.info == "p2") {
    id <- which(ploidy.p1 == dosage.p1[mrk.id])
    g[id, ] <- g[id, ] - ploidy.p1/2
    for(i in phase.conf){
      cat("   Conf.", i,":")
      w <- calc_haploprob_biallelic_single(PH = ph[[i]]$p2,
                                           G = g,
                                           rf = ph[[i]]$rf,
                                           err = ph[[i]]$error)
      ph[[i]]$haploprob <- cbind(rep(2, ploidy.p2*n.ind), w)
    }
    cat("Done with haplotype probabilities\n")
    return(ph)
  }
  else {stop("it should not get here")}
}
