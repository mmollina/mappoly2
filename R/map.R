#' Multi-Locus Map Estimation
#'
#' Performs multi-locus map estimation on mappoly2 data objects. This function supports
#' both sequential and parallel processing for handling large datasets efficiently.
#'
#' @param x An object representing mappoly2 data, typically containing genetic information
#'          and maps for linkage groups.
#' @param lg Optional; a vector of linkage group indices to be processed. If NULL, all
#'           linkage groups in the data object are considered.
#' @param type Character vector indicating the type of mapping to perform. Options include
#'             "mds" (multi-dimensional scaling) and "genome". Default is c("mds", "genome").
#' @param phase.conf A configuration parameter for phase determination.
#' @param rf Recombination fraction, used in the mapping calculations.
#' @param error Error tolerance in the mapping process.
#' @param verbose Logical; if TRUE, detailed progress information will be printed during
#'                processing. Default is TRUE.
#' @param tol Tolerance level for the mapping algorithm.
#' @param ret_H0 Logical; if TRUE, some hypothesis testing result is returned.
#' @param ncpus Integer; specifies the number of cores to use for parallel processing
#'
#' @return Returns an updated mappoly2 data object with mapped linkage groups.
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
#'
#' @importFrom parallel makeCluster stopCluster detectCores
#' @seealso \link[mappoly2]{mapping_one} for the function used in mapping calculations.
mapping <- function(x,
                    lg = NULL,
                    type = c("mds", "genome", "custom"),
                    parent = c("p1p2","p1","p2"),
                    recompute.from.pairwise = FALSE,
                    phase.conf = "all",
                    rf = NULL,
                    error = 0.0,
                    ncpus = 1,
                    verbose = TRUE,
                    tol = 10e-4,
                    ret_H0 = FALSE)
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
  cat("Multi-locus map estimation\n")

  # Assessing phase availability
  # Two-point based
  has.rf.phase <- sapply(x$maps[y$lg], function(x) sapply(x[[y$type]][3:5], function(x) !is.null(x$rf.phase)))
  # HMM screened
  has.hmm.phase <- sapply(x$maps[y$lg], function(x) sapply(x[[y$type]][3:5], function(x) !is.null(x$hmm.phase)))

  # Checking for minimal phase information
  if(all(!has.rf.phase) &  all(!has.hmm.phase))
    stop("Provide a pre-phased sequence.")

  ## Selecting markers based on input
  if(recompute.from.pairwise){
    u <- t(has.rf.phase[parent, , drop = FALSE])
  } else {
    u <- cbind(has.rf.phase[parent,],
               has.hmm.phase[parent,])
  }
  v <- apply(u, 1, any)
  if(all(!v))
    stop(paste("Provide a pre-phased sequence for groups", paste(names(v), collapse = " ")))
  if(any(!v))
    warning(paste("Provide a pre-phased sequence for groups", paste(names(v[!v]), collapse = " ")))
  if(recompute.from.pairwise)
    p <- apply(u[v,,drop= FALSE], 1, function(x) ifelse(x[1], "rf.phase", NA))
  else
    p <- apply(u[v,], 1, function(x) ifelse(x[2], "hmm.phase", "rf.phase"))

  ## Gathering data for parallel processing
  mapData <- vector("list", length(p))
  names(mapData) <- names(p)
  for(i in names(p)){
    mrk.id <- mappoly2:::get_markers_from_phased_sequence(x, i, y$type, parent, phase = p[i])
    mrk.id <- get_info_markers(mrk.id[[i]], x, parent) ## Double checking marker information
    gtemp <- g[mrk.id, ind.names]
    ph <- x$maps[[i]][[y$type]][[parent]][[p[i]]]
    mapData[[i]] <- list(gtemp = gtemp,
                         ph = ph,
                         ploidy.p1 = ploidy.p1,
                         ploidy.p2 = ploidy.p2,
                         dosage.p1= dosage.p1[mrk.id],
                         dosage.p2 = dosage.p2[mrk.id],
                         info.parent = parent,
                         phase.conf = phase.conf,
                         rf = rf,
                         error = error,
                         verbose = verbose,
                         tol = tol,
                         ret_H0 = ret_H0)
  }
  mapFunc <- function(data) {
    if(data$verbose) cat("Processing linkage group\n")
    return(mappoly2:::mapping_one(g = data$gtemp,
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
  if(ncpus > 1) {
    cl <- makeCluster(ncpus)
    results <- parLapply(cl, mapData, mapFunc)
    stopCluster(cl)
  } else {
    results <- lapply(mapData, mapFunc)
  }
  for(i in names(p))
    x$maps[[i]][[y$type]][[parent]][["hmm.phase"]] <- results[[i]]
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
      w <- mappoly2:::est_hmm_map_biallelic(PH = list(ph[[i]]$p1[mrk.id, ],
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
      w <- mappoly2:::est_hmm_map_biallelic_single(PH = ph[[i]]$p1[mrk.id, ],
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
    }
    cat("Done with map estimation\n")
    return(sort_phase(ph))
  }
  else {stop("it should not get here")}
}

#' Efficiently phases unprocessed markers in a given phased
#' genetic map, circumventing the need for complete
#' HMM-based map recomputation.
#'
#' @param void internal function
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
augment_phased_map <- function(x,
                               lg = NULL,
                               type = c("mds", "genome"),
                               thresh.LOD.ph = 5,
                               thresh.LOD.rf = 5,
                               thresh.rf = 0.5,
                               max.phases = 5,
                               thresh.LOD.ph.to.insert = 10,
                               thresh.rf.to.insert = NULL,
                               tol = 10e-4,
                               verbose = TRUE){
  # Extract the linkage group and type information from the input object
  y <- mappoly2:::parse_lg_and_type(x, lg, type)

  # Get all markers for the specified linkage group and type
  mrk.all.lg <- mappoly2:::get_markers_from_ordered_sequence(x, y$lg, y$type, "p1p2")

  # Generate a list of phased map information for each linkage group
  p1p2.map <- lapply(x$maps[y$lg], function(map_item) map_item[[y$type]]$p1p2$hmm.phase[[1]])

  # Check if haploprob is available in the phased map
  assert_that(!is.null(p1p2.map[[i]]$haploprob),
              msg = "Compute haplotype probability for initial sequence")

  # Extract individual names from the screened data
  ind.names <- x$data$screened.data$ind.names
  n.ind <- length(ind.names)  # Number of individuals

  # Extract genotype dosage information and replace NA values with -1
  g <- x$data$geno.dose
  g[is.na(g)] <- -1  # Handling missing data in genotype dosages

  # Get ploidy and dosage information for both parent 1 and parent 2
  ploidy.p1 <- x$data$ploidy.p1
  ploidy.p2 <- x$data$ploidy.p2
  dosage.p1 <- x$data$dosage.p1
  dosage.p2 <- x$data$dosage.p2

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

  for(i in names(mrk.all.lg)){

    # Extract positioned and unpositioned markers
    mrk.pos <- rownames(p1p2.map[[i]]$p1)  # Positioned markers
    mrk.id <- setdiff(mrk.all.lg[[i]], mrk.pos)  # Markers to be positioned

    # Check if there are any markers to be positioned
    if(length(mrk.id) == 0) {
      return(p1p2.map[[i]])  # No unpositioned markers found
    }

    # Function to perform two-point phasing for a parent
    perform_phasing <- function(dosage, map, M, parent) {
      dose.vec <- dosage[mrk.id]
      InitPh <- map
      S <- M[[paste0("Sh.", parent)]][mrk.id, mrk.pos]
      mappoly2:::phasing_one(mrk.id, dose.vec, S, InitPh, verbose = FALSE)
    }

    # Two-point phasing for both parents
    L1 <- perform_phasing(dosage.p1, p1p2.map[[i]]$p1, M[[i]], "p1")
    L2 <- perform_phasing(dosage.p2, p1p2.map[[i]]$p2, M[[i]], "p2")

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
      stop("No markers were selected for 'max.phases' = ", max.phases,
           "\n'max.phases' should be at least ", min(n.conf))
    }

    # Update phase information for selected markers
    L1 <- L1[mrk.sel]
    L2 <- L2[mrk.sel]
    mrk.id <- mrk.id[mrk.sel]

    # Create pedigree matrix
    pedigree <- matrix(rep(c(1, 2, ploidy.p1, ploidy.p2, 1), n.ind),
                       nrow = n.ind, byrow = TRUE)

    # Find flanking markers
    flanking <- mappoly2:::find_flanking_markers(mrk.all.lg[[i]], mrk.pos, mrk.id)
    phasing_results <- vector("list", length(flanking))
    names(phasing_results) <- names(flanking)

    # Initialize progress bar if verbose mode is enabled
    if(verbose) pb <- txtProgressBar(min = 0, max = length(L1), style = 3)

    # Iterating over each set of phasing results
    for(j in 1:length(L1)) {
      # Extract genotype data for the current marker
      G <- g[[i]][mrk.id[j], , drop = TRUE]
      u <- match(unlist(flanking[[mrk.id[j]]]), mrk.pos)

      # Determine the position of the marker and set homolog probabilities
      if(is.na(u)[1]) {  # Marker at the beginning of the linkage group
        homolog_prob <- as.matrix(p1p2.map[[i]]$haploprob[, c(na.omit(u), na.omit(u) + 1) + 3])
        idx <- c(1, 0, 2)
      } else if(is.na(u)[2]) {  # Marker at the end of the linkage group
        homolog_prob <- as.matrix(p1p2.map[[i]]$haploprob[, c(na.omit(u) - 1, na.omit(u)) + 3])
        idx <- c(0, 2, 1)
      } else {  # Marker in the middle of the linkage group
        homolog_prob <- as.matrix(p1p2.map[[i]]$haploprob[, u + 3])
        idx <- c(0, 1, 2)
      }

      # Initialize variables for phasing computations
      w2 <- w1 <- NULL
      z <- vector("list", nrow(L1[[j]]) * nrow(L2[[j]]))
      count <- 1

      # Nested loops for computing phasing results
      for(l in 1:nrow(L1[[j]])) {
        for(k in 1:nrow(L2[[j]])) {
          PH <- list(L1[[j]][l, ], L2[[j]][k, ])
          z[[count]] <- mappoly2:::est_hmm_map_biallelic_insert_marker(PH, G, pedigree, homolog_prob,
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
      if(verbose) setTxtProgressBar(pb, j)
    }

    # Close the progress bar if verbose mode is enabled
    if(verbose) close(pb)


    # Selecting the list of phasing results based on loglike criteria
    selected.list <- phasing_results[sapply(phasing_results, function(x) {
      length(x$loglike) == 1 || x$loglike[2] > thresh.LOD.ph.to.insert
    })]

    # Set default value for thresh.rf.to.insert if it is NULL
    if(is.null(thresh.rf.to.insert)) {
      thresh.rf.to.insert <- max(p1p2.map[[i]]$rf)
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
      cur.mrk <- rownames(p1p2.map[[i]]$p1)
      pos <- mappoly2:::find_flanking_markers(mrk.all.lg[[i]], cur.mrk, j)

      # Skip if no flanking markers are found
      if(length(unlist(pos)) == 0) {
        next()
      }

      # Inserting markers at the beginning, end, or middle of the linkage group
      if(is.na(pos[[1]]$preceding)) {  # Beginning
        p1p2.map[[i]]$p1 <- rbind(selected.list[[j]]$phases$p1[1,], p1p2.map[[i]]$p1)
        p1p2.map[[i]]$p2 <- rbind(selected.list[[j]]$phases$p2[1,], p1p2.map[[i]]$p2)
        rownames(p1p2.map[[i]]$p1) <- rownames(p1p2.map[[i]]$p2) <- c(j, cur.mrk)
      } else if(is.na(pos[[1]]$succeeding)) {  # End
        p1p2.map[[i]]$p1 <- rbind(p1p2.map[[i]]$p1, selected.list[[j]]$phases$p1[1,])
        p1p2.map[[i]]$p2 <- rbind(p1p2.map[[i]]$p2, selected.list[[j]]$phases$p2[1,])
        rownames(p1p2.map[[i]]$p1) <- rownames(p1p2.map[[i]]$p2) <- c(cur.mrk, j)
      } else {  # Middle
        preceding <- cur.mrk[1:match(pos[[1]]$preceding, cur.mrk)]
        succeeding <- cur.mrk[(match(pos[[1]]$succeeding, cur.mrk)):length(cur.mrk)]
        p1p2.map[[i]]$p1 <- rbind(p1p2.map[[i]]$p1[preceding,], selected.list[[j]]$phases$p1[1,], p1p2.map[[i]]$p1[succeeding,])
        p1p2.map[[i]]$p2 <- rbind(p1p2.map[[i]]$p2[preceding,], selected.list[[j]]$phases$p2[1,], p1p2.map[[i]]$p2[succeeding,])
        rownames(p1p2.map[[i]]$p1) <- rownames(p1p2.map[[i]]$p2) <- c(preceding, j, succeeding)
      }
    }

    # Update the hmm phase in the original x$maps object
    x$maps[[i]][[y$type]]$p1p2$hmm.phase[[1]] <- p1p2.map[[i]]

  }
  return(x)
}

#'This function merges two genetic maps built from markers that are exclusively
#'informative in isolated parents, facilitating unified analysis and visualization
#'of distinct genetic data.
#'
#' @param void internal function
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
  y <- mappoly2:::parse_lg_and_type(x,lg,type)
  # HMM screened
  has.hmm.phase <- sapply(x$maps[y$lg], function(x) sapply(x[[y$type]][3:4], function(x) !is.null(x$hmm.phase)))
  if(any(!has.hmm.phase))
    stop("Provide an HMM screened sequence.")
  # Gathering P1 and P2 maps
  p1.map <- lapply(x$maps[y$lg], function(x) x[[y$type]]$p1$hmm.phase[[1]])
  p2.map <- lapply(x$maps[y$lg], function(x) x[[y$type]]$p2$hmm.phase[[1]])
  # Gathering marker order
  ord.mrk.id <- mappoly2:::get_markers_from_ordered_sequence(x, y$lg, y$type, "p1p2")

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

