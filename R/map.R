#' Multi-Locus Map Estimation for mappoly2 Data Objects
#'
#' This function performs multi-locus map estimation on mappoly2 data objects, supporting
#' both sequential and parallel processing. It is designed to efficiently handle large
#' datasets, making it suitable for complex genetic mapping tasks.
#'
#' @param x An object of class \code{mappoly2.sequence}, representing the genetic map data.
#' @param lg Optional vector specifying linkage group indices to be processed. If NULL,
#'           all linkage groups in the data object are considered.
#' @param type A character vector indicating the type of mapping to perform. Options
#'             include "mds" (multi-dimensional scaling), "genome", and "custom".
#'             Default is c("mds", "genome", "custom").
#' @param parent A character vector specifying the parent or parents to be considered
#'               in the mapping process. Options are "p1p2" (both parents), "p1" (first parent),
#'               and "p2" (second parent). Default is c("p1p2", "p1", "p2").
#' @param recompute.from.pairwise A logical value indicating whether to recompute the
#'                                map from pairwise data.
#' @param phase.conf A configuration parameter for phase determination.
#' @param rf Recombination fraction, used in the mapping calculations.
#' @param error Error tolerance in the mapping process.
#' @param ncpus Integer specifying the number of CPU cores for parallel processing.
#'              Default is 1 (sequential processing).
#' @param verbose Logical value; if TRUE, detailed progress information is printed during
#'                processing. Default is TRUE.
#' @param tol Tolerance level for the mapping algorithm.
#' @param ret_H0 Logical; if TRUE, hypothesis testing results are returned.
#'
#' @return Returns an updated \code{mappoly2.sequence} data object with mapped linkage groups.
#'
#' @details The function processes each specified linkage group, performing phase determination
#'          and map estimation based on the provided parameters. It can utilize parallel
#'          processing to enhance performance on large datasets.
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @importFrom parallel makeCluster stopCluster detectCores
#' @export
mapping <- function(x, lg = NULL, type = c("mds", "genome", "custom"),
                    parent = c("p1p2","p1","p2"),
                    recompute.from.pairwise = FALSE,
                    phase.conf = "all", rf = NULL,
                    error = 0.0, ncpus = 1, verbose = TRUE,
                    tol = 10e-4, ret_H0 = FALSE) {
  # Parse the linkage group and map type from 'x'
  y <- parse_lg_and_type(x, lg, type)

  # Validate the input sequence and number of CPUs
  assert_that(is.mappoly2.sequence(x))
  assert_that(is.numeric(ncpus))

  # Adjust the number of CPUs to not exceed available cores
  ncpus <- min(ncpus, detectCores())

  # Match the 'parent' argument to its possible values
  parent <- match.arg(parent)

  # Prepare the genetic data for mapping
  g <- x$data$geno.dose
  g[is.na(g)] <- -1
  ploidy.p1 <- x$data$ploidy.p1
  ploidy.p2 <- x$data$ploidy.p2
  dosage.p1 <- x$data$dosage.p1
  dosage.p2 <- x$data$dosage.p2
  ind.names <- x$data$screened.data$ind.names
  cat("Multi-locus map estimation\n")

  # Assess phase availability in the maps
  has.rf.phase <- sapply(x$maps[y$lg], function(x) sapply(x[[y$type]][3:5], function(x) !is.null(x$rf.phase)))
  has.hmm.phase <- sapply(x$maps[y$lg], function(x) sapply(x[[y$type]][3:5], function(x) !is.null(x$hmm.phase)))

  # Check if there is minimal phase information for mapping
  if(all(!has.rf.phase) & all(!has.hmm.phase)) {
    stop("Provide a pre-phased sequence.")
  }

  # Prepare marker selection based on input and available phase information
  if(recompute.from.pairwise){
    u <- t(has.rf.phase[parent, ,drop = FALSE])
  } else {
    u <- cbind(has.rf.phase[parent,], has.hmm.phase[parent,])
    rownames(u) <- colnames(has.hmm.phase)
  }
  v <- apply(u, 1, any)
  if(all(!v)) {
    stop(paste("Provide a pre-phased sequence for groups", paste(names(v), collapse = " ")))
  }
  if(any(!v)) {
    warning(paste("Provide a pre-phased sequence for groups", paste(names(v[!v]), collapse = " ")))
  }
  if(recompute.from.pairwise) {
    p <- apply(u[v,,drop= FALSE], 1, function(x) ifelse(x[1], "rf.phase", NA))
  } else {
    p <- apply(u[v,,drop= FALSE], 1, function(x) ifelse(x[2], "hmm.phase", "rf.phase"))
  }

  # Prepare data for each map for parallel processing
  mapData <- vector("list", length(p))
  names(mapData) <- names(p)
  for(i in names(p)) {
    mrk.id <- get_markers_from_phased_sequence(x, i, y$type, parent, phase = p[i])
    mrk.id <- get_info_markers(mrk.id[[i]], x, parent)  # Double checking marker information
    gtemp <- g[mrk.id, ind.names]
    ph <- x$maps[[i]][[y$type]][[parent]][[p[i]]]
    mapData[[i]] <- list(gtemp = gtemp,
                         ph = ph,
                         ploidy.p1 = ploidy.p1,
                         ploidy.p2 = ploidy.p2,
                         dosage.p1 = dosage.p1[mrk.id],
                         dosage.p2 = dosage.p2[mrk.id],
                         info.parent = parent,
                         phase.conf = phase.conf,
                         rf = rf,
                         error = error,
                         verbose = verbose,
                         tol = tol,
                         ret_H0 = ret_H0)
  }

  # Define the mapping function for each map
  mapFunc <- function(data) {
    if(data$verbose) cat("Processing linkage group\n")
    return(mapping_one(g = data$gtemp,
                       ph = data$ph,
                       ploidy.p1 = data$ploidy.p1,
                       ploidy.p2 = data$ploidy.p2,
                       dosage.p1 = data$dosage.p1,
                       dosage.p2 = data$dosage.p2,
                       info.parent = data$info.parent,
                       phase.conf = data$phase.conf,
                       rf = data$rf,
                       error = data$error,
                       verbose = data$verbose,
                       tol = data$tol,
                       ret_H0 = data$ret_H0))
  }

  # Execute the mapping in parallel
  if(ncpus > 1) {
    os_type <- Sys.info()["sysname"]
    if (os_type == "Windows") {
      # Windows OS: Use parLapply
      cl <- makeCluster(ncpus)
      on.exit(stopCluster(cl))
      clusterExport(cl, varlist = c("mapData", "mapFunc", "mapping_one"), envir = environment())
      results <- parLapply(cl, mapData, mapFunc)
    } else {
      # Non-Windows OS: Use mclapply
      results <- mclapply(mapData, mapFunc, mc.cores = ncpus)
    }
  } else {
    # Single-core execution: Use lapply
    results <- lapply(mapData, mapFunc)
  }

  # Update the maps with the results
  for(i in names(p)) {
    x$maps[[i]][[y$type]][[parent]][["hmm.phase"]] <- results[[i]]
  }

  return(x)
}


mapping_one <- function(g,
                        ph,
                        ploidy.p1,
                        ploidy.p2,
                        dosage.p1,
                        dosage.p2,
                        info.parent = "p1p2",
                        phase.conf = "all",
                        rf = NULL,
                        error = 0.0,
                        verbose = TRUE,
                        tol = 10e-4,
                        ret_H0 = FALSE){
  if(all(phase.conf == "all"))
    phase.conf <- 1:length(ph)
  assert_that(all(phase.conf%in%1:length(ph)),
              msg = "invalid phases specified in 'phase.conf'")
  n.ind <- ncol(g)
  mrk.id <- rownames(g)
  if(is.null(rf))
    rf <- rep(0.01, nrow(g) - 1)
  assert_that(length(rf) == nrow(g) - 1)
  if (info.parent == "p1p2"){ ###FIXME: include detect_parent_info
    for(i in phase.conf){
      cat("   Conf.", i,":")
      pedigree <- matrix(rep(c(1,2,ploidy.p1,ploidy.p2, 1),n.ind),
                         nrow = n.ind,
                         byrow = TRUE)
      w <- est_hmm_map_biallelic(PH = list(ph[[i]]$p1[mrk.id, ],
                                           ph[[i]]$p2[mrk.id, ]),
                                 G = g,
                                 pedigree = pedigree,
                                 rf = rf,
                                 err = error,
                                 verbose = verbose,
                                 detailed_verbose = FALSE,
                                 tol = tol,
                                 ret_H0 = ret_H0)
      ph[[i]]$loglike <- w[[1]]
      ph[[i]]$rf <- w[[2]]
      ph[[i]]$error <- error
      ph[[i]]$p1 <- ph[[i]]$p1[mrk.id, ]
      ph[[i]]$p2 <- ph[[i]]$p2[mrk.id, ]
    }
    cat("Done with map estimation\n")
    #return(ph)
    return(sort_phase(ph))
  }
  else if(info.parent == "p1"){
    id <- which(ploidy.p2 == dosage.p2[mrk.id])
    g[id, ] <- g[id, ] - ploidy.p2/2
    for(i in phase.conf){
      cat("   Conf.", i,":")
      w <- est_hmm_map_biallelic_single(PH = ph[[i]]$p1[mrk.id, ],
                                        G = g,
                                        rf = rf,
                                        err = error,
                                        verbose = verbose,
                                        detailed_verbose = FALSE,
                                        tol = tol,
                                        ret_H0 = ret_H0)
      ph[[i]]$loglike <- w[[1]]
      ph[[i]]$rf <- w[[2]]
      ph[[i]]$error <- error
      ph[[i]]$p1 <- ph[[i]]$p1[mrk.id, ]
      ph[[i]]$p2 <- ph[[i]]$p2[mrk.id, ]
    }
    cat("Done with map estimation\n")
    return(sort_phase(ph))
  }
  else if(info.parent == "p2") {
    id <- which(ploidy.p1 == dosage.p1[mrk.id])
    g[id, ] <- g[id, ] - ploidy.p1/2
    for(i in phase.conf){
      cat("   Conf.", i,":")
      w <- est_hmm_map_biallelic_single(PH = ph[[i]]$p2[mrk.id, ],
                                        G = g,
                                        rf = rf,
                                        err = error,
                                        verbose = verbose,
                                        detailed_verbose = FALSE,
                                        tol = tol,
                                        ret_H0 = ret_H0)
      ph[[i]]$loglike <- w[[1]]
      ph[[i]]$rf <- w[[2]]
      ph[[i]]$error <- error
      ph[[i]]$p1 <- ph[[i]]$p1[mrk.id, ]
      ph[[i]]$p2 <- ph[[i]]$p2[mrk.id, ]
    }
    cat("Done with map estimation\n")
    return(sort_phase(ph))
  }
  else {stop("it should not get here")}
}

#' Augment a Phased Genetic Map with Unprocessed Markers
#'
#' This function efficiently phases unprocessed markers in a given phased genetic
#' map, circumventing the need for complete HMM-based recomputation of the map.
#' It is designed to update a genetic map with new marker data while maintaining
#' the integrity and structure of the existing map.
#'
#' @param x An object of class \code{mappoly2.sequence}
#' @param lg Optional vector specifying the linkage groups to be processed.
#'           If NULL, all linkage groups in the object are considered.
#' @param type The type of genetic data to be processed, either 'mds' or 'genome'.
#' @param ncpus The number of CPU cores to use for parallel processing.
#' @param thresh.LOD.ph Threshold for the LOD (Logarithm of the Odds) score for
#'                      phasing.
#' @param thresh.LOD.rf Threshold for the LOD score for recombination fractions.
#' @param thresh.rf Threshold for recombination fraction.
#' @param max.phases The maximum number of phase configurations allowed for a marker.
#' @param thresh.LOD.ph.to.insert Threshold LOD score for inserting a marker into the map.
#' @param thresh.dist.to.insert Optional threshold for marker distance (in cM) when
#'                            inserting markers. If NULL, the maximum distance in the map is used.
#' @param reestimate.hmm Logical flag indicating whether to reestimate the HMM
#'                       (Hidden Markov Model) after marker insertion.
#' @param tol Tolerance level for numerical computations.
#' @param final.tol Final tolerance level for HMM reestimation.
#' @param final.error Final error level for HMM reestimation.
#' @param verbose Logical flag indicating whether to print detailed output during
#'                function execution.
#'
#' @return Returns the updated genetic map object with newly phased markers.
#'
#' @details The function works by first identifying unprocessed markers in the
#'          genetic map and then using a phasing algorithm to integrate these
#'          markers into the existing map. It can handle large datasets and is
#'          optimized for performance with options for parallel processing.
#'
#' @importFrom parallel makeCluster stopCluster
#' @importFrom assertthat assert_that
#' @export
augment_phased_map <- function(x,
                               lg = NULL,
                               type = c("mds", "genome"),
                               ncpus = 1,
                               thresh.LOD.ph = 5,
                               thresh.LOD.rf = 5,
                               thresh.rf = 0.5,
                               max.phases = 5,
                               thresh.LOD.ph.to.insert = 10,
                               thresh.dist.to.insert = NULL,
                               reestimate.hmm = TRUE,
                               tol = 10e-3,
                               final.tol = 10e-4,
                               final.error = NULL,
                               verbose = TRUE){
  # Extract the linkage group and type information from the input object
  y <- parse_lg_and_type(x, lg, type)

  # Get all markers for the specified linkage group and type
  mrk.all.lg <- get_markers_from_ordered_sequence(x, y$lg, y$type, "p1p2")

  # Generate a list of phased map information for each linkage group
  p1p2.map <- lapply(x$maps[y$lg], function(map_item) map_item[[y$type]]$p1p2$hmm.phase[[1]])

  # Check if haploprob is available in the phased map
  assert_that(all(sapply(p1p2.map, function(x) !is.null(x$haploprob))),
              msg = "Compute haplotype probability for initial sequence")

  # Extract individual names from the screened data
  ind.names <- x$data$screened.data$ind.names
  n.ind <- length(ind.names)  # Number of individuals

  if(is.null(final.error))
    final.error <- x$maps[[1]][[y$type]]$p1p2$hmm.phase[[1]]$error

  # Extract genotype dosage information and replace NA values with -1
  g <- x$data$geno.dose[,ind.names]
  g[is.na(g)] <- -1  # Handling missing data in genotype dosages

  # Get ploidy and dosage information for both parent 1 and parent 2
  ploidy.p1 <- x$data$ploidy.p1
  ploidy.p2 <- x$data$ploidy.p2
  dosage.p1 <- x$data$dosage.p1
  dosage.p2 <- x$data$dosage.p2

  thresh.rf.to.insert <- mf_h(thresh.dist.to.insert)

  # Retrieving how many alternate alleles share homologs based
  # on pairwise linkage analysis
  M <- lapply(mrk.all.lg,
              function(mrk.seq) filter_rf_matrix(x$data,
                                                 type = "sh",
                                                 thresh.LOD.ph,
                                                 thresh.LOD.rf,
                                                 thresh.rf,
                                                 mrk.names = mrk.seq))
  # splitting data into linkage groups
  g <- lapply(mrk.all.lg, function(x) g[x, ])

  # Preparing data to serial or parallel submission
  mapData <- vector("list", length(mrk.all.lg))
  names(mapData) <- names(mrk.all.lg)

  for(i in names(mrk.all.lg)){
    mapData[[i]] <- list(map = p1p2.map[[i]],
                         mrk = mrk.all.lg[[i]],
                         mat = M[[i]],
                         geno = g[[i]],
                         max.phases = max.phases,
                         ploidy.p1 = ploidy.p1,
                         ploidy.p2 = ploidy.p2,
                         dosage.p1 = dosage.p1,
                         dosage.p2 = dosage.p2,
                         tol = tol,
                         thresh.LOD.ph.to.insert = thresh.LOD.ph.to.insert,
                         thresh.rf.to.insert = thresh.rf.to.insert,
                         verbose = verbose,
                         n.ind = n.ind)
  }
  if(ncpus > 1) {
    cl <- makeCluster(ncpus)
    mapResult <- parLapply(cl, mapData, function(x) augment_phased_map_one(x$map, x$mrk, x$mat,
                                                                           x$geno,
                                                                           x$max.phases,
                                                                           x$ploidy.p1,
                                                                           x$ploidy.p2,
                                                                           x$dosage.p1,
                                                                           x$dosage.p2,
                                                                           x$tol,
                                                                           x$thresh.LOD.ph.to.insert,
                                                                           x$thresh.rf.to.insert,
                                                                           x$verbose,
                                                                           x$n.ind))
    stopCluster(cl)
  } else {
    mapResult <- lapply(mapData, function(x) augment_phased_map_one(x$map, x$mrk, x$mat,
                                                                    x$geno,
                                                                    x$max.phases,
                                                                    x$ploidy.p1,
                                                                    x$ploidy.p2,
                                                                    x$dosage.p1,
                                                                    x$dosage.p2,
                                                                    x$tol,
                                                                    x$thresh.LOD.ph.to.insert,
                                                                    x$thresh.rf.to.insert,
                                                                    x$verbose,
                                                                    x$n.ind))
  }
  for(i in names(mrk.all.lg)){
    # Update the hmm phase in the original x$maps object
    x$maps[[i]][[y$type]]$p1p2$hmm.phase[[1]] <- mapResult[[i]]

  }
  if(reestimate.hmm){
    cat("\nReestimating multilocus map ...\n")
    x <- mapping(x,
                 lg = y$lg,
                 parent = "p1p2",
                 type = y$type,
                 tol = final.tol,
                 ncpus = ncpus,
                 error = final.error,
                 verbose = FALSE)
  }
  return(x)
}


augment_phased_map_one <- function(map, mrk, mat, geno, max.phases,
                                   ploidy.p1, ploidy.p2, dosage.p1,
                                   dosage.p2, tol, thresh.LOD.ph.to.insert,
                                   thresh.rf.to.insert, verbose, n.ind){
  # Extract positioned and unpositioned markers
  mrk.pos <- rownames(map$p1)  # Positioned markers
  mrk.id <- setdiff(mrk, mrk.pos)  # Markers to be positioned

  # Check if there are any markers to be positioned
  if(length(mrk.id) == 0) {
    return(map)  # No unpositioned markers found
  }

  # Function to perform two-point phasing for a parent
  perform_phasing <- function(dosage, map, M, parent, mrk.id, mrk.pos) {
    dose.vec <- dosage[mrk.id]
    InitPh <- map
    S <- M[[paste0("Sh.", parent)]][mrk.id, mrk.pos, drop = FALSE]
    phasing_one(mrk.id, dose.vec, S, InitPh, verbose = FALSE)
  }

  # Two-point phasing for both parents
  L1 <- perform_phasing(dosage.p1, map$p1, mat, "p1", mrk.id, mrk.pos)
  L2 <- perform_phasing(dosage.p2, map$p2, mat, "p2", mrk.id, mrk.pos)

  # Selecting phase configurations
  n.conf <- sapply(L1, nrow) * sapply(L2, nrow)
  if(verbose) {
    cat("Distribution of phase configurations:\n")
    temp <- as.data.frame(table(n.conf))
    colnames(temp) <- c("n. phase. conf.", "frequency")
    print(temp)
  }

  # Handle cases with no selected markers
  mrk.sel <- which(n.conf <= max.phases)
  if(length(mrk.sel) == 0) {
    warning("No markers were selected for 'max.phases' = ", max.phases,
            "\n increasing 'max.phases' to ", min(n.conf) + 1)
    max.phases <- min(n.conf) + 1
    mrk.sel <- which(n.conf <= max.phases)
  }

  # Update phase information for selected markers
  L1 <- L1[mrk.sel]
  L2 <- L2[mrk.sel]
  mrk.id <- mrk.id[mrk.sel]

  # Create pedigree matrix
  pedigree <- matrix(rep(c(1, 2, ploidy.p1, ploidy.p2, 1), n.ind),
                     nrow = n.ind, byrow = TRUE)

  # Find flanking markers
  flanking <- find_flanking_markers(mrk, mrk.pos, mrk.id)
  phasing_results <- vector("list", length(flanking))
  names(phasing_results) <- names(flanking)

  # Initialize progress bar if verbose mode is enabled
  if(verbose) pb <- utils::txtProgressBar(min = 0, max = length(L1), style = 3)

  # Iterating over each set of phasing results
  for(j in 1:length(L1)) {
    # Extract genotype data for the current marker
    G <- geno[mrk.id[j], , drop = TRUE]
    u <- match(unlist(flanking[[mrk.id[j]]]), mrk.pos)

    # Determine the position of the marker and set homolog probabilities
    if(is.na(u)[1]) {  # Marker at the beginning of the linkage group
      homolog_prob <- as.matrix(map$haploprob[, c(na.omit(u), na.omit(u) + 1) + 3])
      idx <- c(0, 1, 2)
    } else if(is.na(u)[2]) {  # Marker at the end of the linkage group
      homolog_prob <- as.matrix(map$haploprob[, c(na.omit(u) - 1, na.omit(u)) + 3])
      idx <- c(0, 2, 1)
    } else {  # Marker in the middle of the linkage group
      homolog_prob <- as.matrix(map$haploprob[, u + 3])
      idx <- c(1, 0, 2)
    }

    # Initialize variables for phasing computations
    w2 <- w1 <- NULL
    z <- vector("list", nrow(L1[[j]]) * nrow(L2[[j]]))
    count <- 1

    # Nested loops for computing phasing results
    for(l in 1:nrow(L1[[j]])) {
      for(k in 1:nrow(L2[[j]])) {
        PH <- list(L1[[j]][l, ], L2[[j]][k, ])
        z[[count]] <- est_hmm_map_biallelic_insert_marker(PH, G, pedigree, homolog_prob,
                                                          rf = c(0.01, 0.01), idx, verbose = FALSE,
                                                          detailed_verbose = FALSE, tol = tol, ret_H0 = FALSE)
        w1 <- rbind(w1, L1[[j]][l, ])
        w2 <- rbind(w2, L2[[j]][k, ])
        count <- count + 1
      }
    }

    # Process the phasing results
    v <- sapply(z, function(x) x[[1]])
    v <- max(v) - v
    id <- order(v)
    phasing_results[[mrk.id[j]]] <- list(loglike = v[id],
                                         rf.vec = t(sapply(z[id], function(x) x[[2]])),
                                         phases = list(p1 = w1[id, , drop = FALSE],
                                                       p2 = w2[id, , drop = FALSE]))

    # Update progress bar if verbose mode is enabled
    if(verbose) utils::setTxtProgressBar(pb, j)
  }

  # Close the progress bar if verbose mode is enabled
  if(verbose) close(pb)


  # Selecting the list of phasing results based on loglike criteria
  selected.list <- phasing_results[sapply(phasing_results, function(x) {
    length(x$loglike) == 1 || x$loglike[2] > thresh.LOD.ph.to.insert
  })]

  # Set default value for thresh.rf.to.insert if it is NULL
  if(is.null(thresh.rf.to.insert)) {
    thresh.rf.to.insert <- max(map$rf)
  }

  # Validate thresh.rf.to.insert value
  if(thresh.rf.to.insert < 0 || thresh.rf.to.insert >= 0.5) {
    stop("'thresh.rf.to.insert' parameter must be between 0 and 0.5")
  }

  # Further filtering of phasing results based on recombination frequency
  selected.list <- phasing_results[sapply(phasing_results, function(x) {
    max(x$rf.vec[1,]) <= thresh.rf.to.insert
  })]

  # Iterate over selected list to update map information
  for(j in names(selected.list)) {
    cur.mrk <- rownames(map$p1)
    pos <- find_flanking_markers(mrk, cur.mrk, j)

    # Skip if no flanking markers are found
    if(length(unlist(pos)) == 0) {
      next()
    }

    # Inserting markers at the beginning, end, or middle of the linkage group
    if(is.na(pos[[1]]$preceding)) {  # Beginning
      map$p1 <- rbind(selected.list[[j]]$phases$p1[1,], map$p1)
      map$p2 <- rbind(selected.list[[j]]$phases$p2[1,], map$p2)
      rownames(map$p1) <- rownames(map$p2) <- c(j, cur.mrk)
    } else if(is.na(pos[[1]]$succeeding)) {  # End
      map$p1 <- rbind(map$p1, selected.list[[j]]$phases$p1[1,])
      map$p2 <- rbind(map$p2, selected.list[[j]]$phases$p2[1,])
      rownames(map$p1) <- rownames(map$p2) <- c(cur.mrk, j)
    } else {  # Middle
      preceding <- cur.mrk[1:match(pos[[1]]$preceding, cur.mrk)]
      succeeding <- cur.mrk[(match(pos[[1]]$succeeding, cur.mrk)):length(cur.mrk)]
      map$p1 <- rbind(map$p1[preceding,], selected.list[[j]]$phases$p1[1,], map$p1[succeeding,])
      map$p2 <- rbind(map$p2[preceding,], selected.list[[j]]$phases$p2[1,], map$p2[succeeding,])
      rownames(map$p1) <- rownames(map$p2) <- c(preceding, j, succeeding)
    }
  }
  map$loglike <- NULL
  map$rf <- NULL
  map$error <- NULL
  map$haploprob <- NULL
  return(map)
}


#' Merge Single Parent Genetic Maps
#'
#' This function merges two genetic maps built from markers that are exclusively informative
#' in isolated parents. It facilitates unified analysis and visualization of distinct genetic data.
#'
#' @param x An object of class `mappoly2.sequence` containing maps constructed for each parent separately.
#' @param lg Optional vector specifying the linkage groups to be processed.
#'           If NULL, all linkage groups in `x` are considered.
#' @param type The type of genetic maps to be merged, options include "mds", "genome", or "custom".
#' @param hmm.reconstruction Logical; if TRUE, HMM-based reconstruction of the merged map is performed.
#' @param rf Recombination fraction, used in the merging calculations.
#' @param error Error tolerance in the merging process.
#' @param ncpus Integer specifying the number of CPU cores for parallel processing.
#' @param verbose Logical; if TRUE, progress messages will be printed.
#' @param tol Tolerance level for the merging algorithm.
#' @param ret_H0 Logical; if TRUE, hypothesis testing results are returned.
#'
#' @return Returns the `mappoly2.sequence` object with the merged genetic maps for the specified linkage groups.
#'
#' @details The function merges separate genetic maps for individual parents into a single map
#'          for each linkage group. It handles the alignment and integration of markers from
#'          both parents and optionally performs HMM-based reconstruction of the merged maps.
#'
#' @importFrom assertthat assert_that
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
merge_single_parent_maps <- function(x,
                                     lg = NULL,
                                     type = c("mds", "genome", "custom"),
                                     hmm.reconstruction = TRUE,
                                     rf = NULL,
                                     error = 0.0,
                                     ncpus = 1,
                                     verbose = TRUE,
                                     tol = 10e-4,
                                     ret_H0 = FALSE)
{
  y <- parse_lg_and_type(x,lg,type)
  # HMM screened
  has.hmm.phase <- sapply(x$maps[y$lg], function(x) sapply(x[[y$type]][3:4], function(x) !is.null(x$hmm.phase)))
  if(any(!has.hmm.phase))
    stop("Provide an HMM screened sequence.")
  # Gathering P1 and P2 maps
  p1.map <- lapply(x$maps[y$lg], function(x) x[[y$type]]$p1$hmm.phase[[1]])
  p2.map <- lapply(x$maps[y$lg], function(x) x[[y$type]]$p2$hmm.phase[[1]])
  # Gathering marker order
  ord.mrk.id <- get_markers_from_ordered_sequence(x, y$lg, y$type, "p1p2")

  for(i in names(ord.mrk.id)){
    temp.mrk.id <- character(length(ord.mrk.id[[i]]))
    temp.mrk.id[match(rownames(p1.map[[i]]$p1), ord.mrk.id[[i]])] <- rownames(p1.map[[i]]$p1)
    temp.mrk.id[match(rownames(p2.map[[i]]$p2), ord.mrk.id[[i]])] <- rownames(p2.map[[i]]$p2)
    temp.mrk.id <- temp.mrk.id[temp.mrk.id != ""]
    ph.p1 <- matrix(0, length(temp.mrk.id), x$data$ploidy.p1, dimnames = list(temp.mrk.id, NULL))
    ph.p2 <- matrix(0, length(temp.mrk.id), x$data$ploidy.p2, dimnames = list(temp.mrk.id, NULL))
    ph.p1[rownames(p1.map[[i]]$p1), ] <- p1.map[[i]]$p1
    ph.p2[rownames(p1.map[[i]]$p2), ] <- p1.map[[i]]$p2
    ph.p1[rownames(p2.map[[i]]$p1), ] <- p2.map[[i]]$p1
    ph.p2[rownames(p2.map[[i]]$p2), ] <- p2.map[[i]]$p2
    x$maps[[i]][[type]]$p1p2$hmm.phase[[1]] <- list(p1 = ph.p1,
                                                    p2 = ph.p2,
                                                    loglike = NULL,
                                                    rf = NULL,
                                                    error = NULL,
                                                    haploprob = NULL)
  }
  if(hmm.reconstruction)
    x <- mapping(x, lg = lg, type = type, parent = "p1p2",
                 rf = rf, error = error, ncpus = ncpus,
                 verbose = verbose, tol = tol, ret_H0 = ret_H0)
  return(x)
}



#' Compare Marker Orders from MDS and Genome Maps
#'
#' This function compares two pre-built maps, one using MDS (multidimensional scaling)
#' and the other using genome order, based on their HMM (hidden Markov model) multilocus-based likelihoods.
#' It serves as an objective function to assist in selecting the best map by ensuring
#' that only markers present in both maps are considered for a proper comparison.
#'
#' @param x A `mappoly2.sequence` object containing the pre-built maps.
#' @param parent The parent phase to consider ('p1p2', 'p1', or 'p2').
#' @param ncpus Integer specifying the number of CPU cores for parallel processing.
#'              Default is 1 (sequential processing).
#' @param error Optional error parameter for comparison adjustments.
#' @param verbose Logical; if `TRUE`, prints messages during the function execution.
#' @param tol Tolerance level for likelihood comparison.
#'
#' @return A matrix comparing the HMM likelihoods of markers across MDS and genome maps.
#' @export
compare_order <- function(x,
                          parent = c("p1p2","p1","p2"),
                          ncpus = 1,
                          error = 0.0,
                          verbose = TRUE,
                          tol = 10e-3) {
  # Set parent type
  parent <- match.arg(parent)

  assert_that(is.numeric(ncpus))

  if(ncpus > detectCores())
    ncpus <- detectCores() - 1

  # Handle missing genotype data
  g <- x$data$geno.dose
  g[is.na(g)] <- -1

  # Extract necessary data from the input object
  ploidy.p1 <- x$data$ploidy.p1
  ploidy.p2 <- x$data$ploidy.p2
  dosage.p1 <- x$data$dosage.p1
  dosage.p2 <- x$data$dosage.p2
  ind.names <- x$data$screened.data$ind.names

  # Initialize list to store comparison results for each linkage group
  d <- w <- vector("list", length(x$maps))
  names(d) <- names(w) <- names(x$maps)

  # Checking existence of maps
  has.mds.map <- sapply(x$maps, function(x) !is.null(x[["mds"]][[parent]]$hmm.phase[[1]]$loglike))
  has.genome.map <- sapply(x$maps, function(x) !is.null(x[["genome"]][[parent]]$hmm.phase[[1]]$loglike))
  id <- has.mds.map & has.genome.map

  if(all(!id))
    stop("there are no maps to compare")

  # Loop over each linkage group to compare orders
  for(i in names(w)[id]) {
    # Extract phase and marker data for MDS and genome methods
    ph.mds <- x$maps[[i]][["mds"]][[parent]]$hmm.phase
    mds.id <- rownames(ph.mds[[1]]$p1)
    d.mds <- c(0, cumsum(imf_h(ph.mds[[1]]$rf)))
    names(d.mds) <- mds.id

    ph.gen <- x$maps[[i]][["genome"]][[parent]]$hmm.phase
    gen.id <- rownames(ph.gen[[1]]$p1)
    d.gen <- c(0, cumsum(imf_h(ph.gen[[1]]$rf)))
    names(d.gen) <- gen.id

    # Intersect markers from both methods for comparison
    mds.id2 <- intersect(mds.id, gen.id)
    gen.id2 <- intersect(gen.id, mds.id)

    # Store comparison data in the list
    w[[i]] <- list(ph.mds = ph.mds,
                   ph.gen = ph.gen,
                   rf.mds = mf_h(diff(d.mds[mds.id2])),
                   rf.gen = mf_h(diff(d.gen[gen.id2])),
                   g.mds = g[mds.id2, ind.names],
                   g.gen = g[gen.id2, ind.names])
    d[[i]] <- list(mds = d.mds, genome = d.gen)
  }
  if(ncpus > 1){
    cl <- makeCluster(ncpus)
    # Apply the comparison function in parallel
    result <- parLapply(cl, w, function(x, parent, ploidy.p1, ploidy.p2, dosage.p1,
                                        dosage.p2, error, verbose, tol) {
      compare_order_one_lg(x$g.mds, x$g.gen, x$ph.mds, x$ph.gen, x$rf.mds,
                           x$rf.gen, parent, ploidy.p1, ploidy.p2, dosage.p1,
                           dosage.p2, error = error, verbose = verbose, tol = tol)
    }, parent, ploidy.p1, ploidy.p2, dosage.p1, dosage.p2, error, verbose, tol)

    # Stop the cluster after computation is done
    stopCluster(cl)
  }
  else {
    # Apply the comparison function to each linkage group
    result <- lapply(w, function(x, parent, ploidy.p1, ploidy.p2, dosage.p1,
                                 dosage.p2, error, verbose, tol) {
      compare_order_one_lg(x$g.mds, x$g.gen, x$ph.mds, x$ph.gen, x$rf.mds,
                           x$rf.gen, parent, ploidy.p1, ploidy.p2, dosage.p1,
                           dosage.p2, error = error, verbose = verbose, tol = tol)
    }, parent, ploidy.p1, ploidy.p2, dosage.p1, dosage.p2, error, verbose, tol)
  }
  # Combine results into a matrix
  comp.mat <- NULL
  for(i in names(result)) {
    comp.mat <- rbind(comp.mat, cbind(LG = i, result[[i]]))
  }
  structure(list(comp.mat = comp.mat, maps = d), class = "mappoly2.order.comparison")
}

#' Compare Marker Orders for a Single Linkage Group
#'
#' This internal function compares the marker orders for a single linkage group
#' using two different map construction methods: multidimensional scaling (MDS)
#' and genome order. It calculates the mapping for each method and compares
#' their log likelihoods to determine the best map order.
#'
#' @param g.mds Genotype data for MDS.
#' @param g.gen Genotype data for genome method.
#' @param ph.mds Phase information for MDS.
#' @param ph.gen Phase information for genome method.
#' @param rf.mds Recombination fractions for MDS.
#' @param rf.gen Recombination fractions for genome method.
#' @param parent Parent phase to consider.
#' @param ploidy.p1 Ploidy of parent 1.
#' @param ploidy.p2 Ploidy of parent 2.
#' @param dosage.p1 Dosage information for parent 1.
#' @param dosage.p2 Dosage information for parent 2.å
#' @param error error parameter.
#' @param verbose If `TRUE`, prints messages.
#' @param tol Tolerance level.
#'
#' @return A matrix with comparison results, including number of markers, map lengths,
#' maximum gaps, average distances, and log likelihoods for both methods.
#' @keywords internal
compare_order_one_lg <- function(g.mds, g.gen, ph.mds, ph.gen, rf.mds,
                                 rf.gen, parent, ploidy.p1, ploidy.p2,
                                 dosage.p1, dosage.p2, error = 0.0,
                                 verbose = TRUE, tol = 10e-3){
  mds.map <- mapping_one(g.mds,
                         ph.mds,
                         ploidy.p1,
                         ploidy.p2,
                         dosage.p1,
                         dosage.p2,
                         parent,
                         phase.conf = "all",
                         rf = rf.mds,
                         error = error,
                         verbose = verbose,
                         tol = tol,
                         ret_H0 = TRUE)
  gen.map <- mapping_one(g.gen,
                         ph.gen,
                         ploidy.p1,
                         ploidy.p2,
                         dosage.p1,
                         dosage.p2,
                         parent,
                         phase.conf = "all",
                         rf = rf.gen,
                         error = error,
                         verbose = verbose,
                         tol = tol,
                         ret_H0 = TRUE)
  comp.mat <- sapply(list(mds.map, gen.map), function(x){
    c(nrow(x[[1]]$p1), sum(imf_h(x[[1]]$rf)), max(imf_h(x[[1]]$rf)),
      mean(imf_h(x[[1]]$rf)), x[[1]]$loglike)})
  comp.mat <- t(round(comp.mat, 2))
  l <- list(c("*", ""), c("", "*"))
  comp.mat <- cbind(comp.mat, l[[which.max(comp.mat[,5])]])
  dimnames(comp.mat) <- list(c("mds", "genome"), c("n.mrk", "map length", "max_gap", "ave_dist",  "loglike", ""))
  return(comp.mat)
}



