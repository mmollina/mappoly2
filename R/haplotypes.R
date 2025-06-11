#' Calculate Haplotype Probabilities Using Hidden Markov Models
#'
#' This function calculates haplotypes for each linkage group in a genetic mapping dataset.
#' It utilizes Hidden Markov Models and supports both sequential and parallel processing
#' for efficient computation.
#'
#' @param x An object representing genetic mapping data, typically of a specific class
#'          that stores genetic information.
#' @param lg Optional vector specifying the linkage group indices to process.
#'           If NULL, all linkage groups in `x` are processed.
#' @param type Character vector indicating the type of map to process, either "mds"
#'             or "genome".
#' @param parent Character vector specifying the parent or parents to be considered
#'               in the haplotype calculation. Options are "p1p2" (both parents),
#'               "p1" (first parent), and "p2" (second parent).
#' @param phase.conf A configuration parameter for phase calculation.
#' @param verbose Logical value; if TRUE, progress messages will be printed.
#' @param ncpus The number of CPU cores to use for parallel processing. Defaults to 1.
#'
#' @return The input object `x` with haplotypes calculated for the specified linkage groups.
#'
#' @details The function processes the genetic data to calculate haplotypes for each
#'          specified linkage group using Hidden Markov Models. It can handle large
#'          datasets efficiently by using parallel processing capabilities.
#'
#' @importFrom parallel detectCores makeCluster parLapply stopCluster
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
calc_haplotypes <- function(x, lg = NULL, type = c("mds", "genome"),
                            parent = c("p1p2","p1","p2"), phase.conf = "all",
                            verbose = TRUE, ncpus = 1) {
  # Parse the linkage group and map type
  # if(is.null(lg)){lg <- 1:length(unique(x$data$chrom))}
  if(length(type) > 1){stop("Please select only one type of analysis at the time.", call. = FALSE)}
  y <- parse_lg_and_type(x, lg, type)

  # Match the 'parent' argument to its possible values
  parent <- match.arg(parent)

  # Prepare the genetic data for haplotype calculation
  g <- x$data$geno.dose
  g[is.na(g)] <- -1  # Replace NA with -1
  ploidy.p1 <- x$data$ploidy.p1
  ploidy.p2 <- x$data$ploidy.p2
  dosage.p1 <- x$data$dosage.p1
  dosage.p2 <- x$data$dosage.p2
  ind.names <- x$data$screened.data$ind.names
  
  # Assess the availability of multi-point map data
  has.hmm.map <- sapply(x$maps[y$lg], function(x) sapply(x[[y$type]][3:5], function(x) !is.null(x$hmm.phase[[1]]$loglike)))

  # Check for the presence of phase information
  if(!all(has.hmm.map)) {
    issues.has.hmm.map <- which(!has.hmm.map, arr.ind = TRUE)
    issues.has.hmm.map <- paste( rownames(has.hmm.map)[issues.has.hmm.map[,1]], colnames(has.hmm.map)[issues.has.hmm.map[,2]] ,sep="_")
    stop(paste("Provide an hmm estimated map for", paste(issues.has.hmm.map, collapse = ", ")), call. = FALSE)
  }

  # Select markers based on the availability of phase information
  u <- t(has.hmm.map[parent, , drop = FALSE])
  v <- apply(u, 1, any)
  if(all(!v)) {
    stop(paste("Provide a pre-phased sequence for groups", paste(names(v), collapse = " ")))
  }
  if(any(!v)) {
    warning(paste("Provide a pre-phased sequence for groups", paste(names(v[!v]), collapse = " ")))
  }
  p <- apply(u[v,,drop = FALSE], 1, function(x) ifelse(x[1], "hmm.phase", NA))

  # Prepare data for each map for parallel processing
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

  # Define the haplotype calculation function
  haplotypeFunc <- function(data) {
    return(calc_haplotypes_one(data$g, data$ph, data$ploidy.p1,
                               data$ploidy.p2, data$dosage.p1,
                               data$dosage.p2, data$phase.conf,
                               parent.info = data$parent.info,
                               verbose = data$verbose))
  }

  # Execute the haplotype calculation in parallel
  if(ncpus > 1) {
    os_type <- Sys.info()["sysname"]
    ncpus <- min(ncpus, detectCores())
    if (os_type == "Windows") {
      cl <- makeCluster(ncpus)
      on.exit(stopCluster(cl))
      clusterExport(cl, varlist = c("haplotypeData", "haplotypeFunc", "calc_haplotypes_one"), envir = environment())
      results <- parLapply(cl, haplotypeData, haplotypeFunc)
    } else {
      results <- mclapply(haplotypeData, haplotypeFunc, mc.cores = ncpus)
    }
  } else {
    # Single-core execution: Use lapply
    results <- lapply(haplotypeData, haplotypeFunc)
  }

  # Update the maps with the calculated haplotypes
  for(i in names(p)) {
    x$maps[[i]][[y$type]][[parent]][["hmm.phase"]] <- results[[i]]
  }

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
  assertthat::assert_that(all(phase.conf%in%1:length(ph)),
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
