#' Match Homologs Across Biparental Genetic Maps
#'
#' This internal function matches shared parents across biparental genetic maps and rearranges
#' the homologs to ensure consistency in joint map estimation. If inconsistencies are found,
#' such as different dosages or phasing, the function removes the marker. It operates on a single
#' chromosome.
#'
#' @param x A list of biparental genetic maps.
#' @param par.ord Data frame or matrix indicating the order of parents.
#' @param lg Linkage group to be analyzed.
#' @param pl Ploidy level.
#' @param par.name Name of the parent.
#' @return A list containing the matched homologs, the hclust object, shared markers, and genome positions.
#' @details This function first identifies shared markers across all maps for a given linkage group
#'          and parent. It then checks for consistency in phase and dosage across these markers.
#'          Inconsistencies lead to the removal of the marker. The function is designed for internal
#'          use and handles complex genetic map structures.
#' @importFrom dplyr mutate arrange
#' @importFrom reshape2 melt acast
#' @importFrom stats dist
#' @noRd
#' @keywords internal
match_homologs <- function(map.list,
                           par.ord,
                           lg,
                           pl,
                           par.name){
  # Number of full-sibs with the analyzed parent
  n <- nrow(par.ord)
  names(map.list) <- sapply(map.list, function(x) paste0(x$data$name.p1, "_x_", x$data$name.p2))

  x <- map.list[rownames(par.ord)]

  # Identify shared markers across all maps for a given linkage group and parent
  idn <- Reduce(intersect, lapply(x, function(x) rownames(x$maps[[lg]]$genome$p1p2$hmm.phase[[1]]$p1)))

  if(length(idn) < 4)
    stop("There is less then 4 markers connecting populations with ", par.name)

  # Collect all markers and their positions from each map
  pos.all <- lapply(x, function(x) data.frame(x$data$mrk.names, x$data$genome.pos))
  pos.all <- Reduce(rbind, pos.all)
  pos.all <- pos.all[!duplicated(pos.all[,1]),]
  dimnames(pos.all) <- list(pos.all[,1], c("mrk", "pos"))
  pos.all <- pos.all[unique(unlist(lapply(x, function(x) rownames(x$maps[[lg]]$genome$p1p2$hmm.phase[[1]]$p1)), use.names = FALSE)),,drop = FALSE]

  # Order markers according to genome positions
  pos <- pos.all[order(pos.all$pos),-1, drop = FALSE]

  # Collect phases for each map and parent
  ph.list <- vector("list", n)
  ph.names <- ph.mat <- NULL
  for(i in 1:n){
    ph <- x[[rownames(par.ord)[i]]]$maps[[lg]]$genome$p1p2$hmm.phase[[1]][[par.ord[i,2]]]
    colnames(ph) <- paste0("H", 1:pl, "_Pop", par.ord[i,1], "_P",  par.ord[i,2])
    ph.list[[i]] <- ph
    ph.mat <- rbind(ph.mat, t(ph[idn,]))
    ph.names <- c(ph.names, paste0("Pop", par.ord[i,1], "_P",  par.ord[i,2]))
  }
  names(ph.list) <- ph.names
  value <- pop <- L1 <- NULL
  # Cluster and rearrange homologs if more than one full-sib family is present
  if(n > 1){
    dd <- as.matrix(dist(ph.mat, method = "binary"))
    for(i in 1:n){
      id <- (((i-1)*pl)+1):(pl*i)
      dd[id,id][] <- 1
    }
    dd <- as.dist(dd)
    hc <- hclust(dd, method = "ward.D2")
    homologs <- cutree(hc, k = pl)
    u <- split(names(homologs), as.factor(homologs))
    names(u) <- paste0("h", 1:pl)
    u <- u %>%
      melt %>%
      mutate(homolog = substr(value, 1,1)) %>%
      mutate(pop = substr(value,4,14)) %>%
      arrange(pop, L1) %>%
      acast(pop ~ L1, value.var = "value")

    # Reorganize homologs for consistency
    for(i in names(ph.list))
      ph.list[[i]] <- ph.list[[i]][,u[i,]]

    # Check for consistency across shared markers
    S <- sapply(ph.list, function(x) apply(x[idn,], 1, paste0, collapse = ""))
    id.ph <- which(apply(S, 1, function(x) length(unique(x)) != 1))
    ph.out <- ph.list[[1]]
    ph.out[names(id.ph),][] <- NA
    for(i in 1:n){
      idtemp <- setdiff(rownames(ph.list[[i]]), rownames(ph.out))
      ph.out <- rbind(ph.out, ph.list[[i]][idtemp,])
    }
    remaining <- setdiff(rownames(pos), rownames(ph.out))
    if(length(remaining) > 0)
      ph.out <- rbind(ph.out, matrix(NA, length(remaining), pl, dimnames = list(remaining, colnames(ph.out))))
  } else {
    hc <- NA
    ph.out <- ph.list[[1]]
    for(i in 1:length(x)){
      idtemp <- setdiff(rownames(x[[i]]$maps[[lg]]$genome$p1p2$hmm.phase[[1]]$p1), rownames(ph.out))
      ph.out <- rbind(ph.out, matrix(NA, length(idtemp), pl, dimnames = list(idtemp, NULL)))
    }
  }

  # Rearrange the output according to genome positions
  ph.out <- ph.out[rownames(pos),]
  colnames(pos) <- "geno.pos"

  # Return the final list of outputs
  list(ph = ph.out, hc = hc, shared.mrks = idn, genome.pos = pos)
}


#' Prepare Genetic Map Data for Integration
#'
#' This function prepares multiple biparental genetic maps for integration into
#' a unified multi-population genetic map. It primarily utilizes an internal
#' function `match_homologs` to match and rearrange homologs across shared
#' parents in different populations, ensuring consistency in the integrated
#' map estimation. Markers with inconsistencies, such as varying dosages or phasing,
#' are removed. The function currently supports only genome-ordered maps
#' and operates on a single chromosome.
#'
#' @param x A list of genetic map objects.
#' @param lg Linkage group to be analyzed, default is NULL.
#' @param type Type of the map, either "genome" or "mds" (default is "genome").
#' @param verbose Logical, indicating whether to show detailed messages (default is TRUE).
#' @return A list of prepared data for building an integrated genetic map.
#' @importFrom tibble column_to_rownames
#' @export
prepare_to_integrate <- function(x,
                                 lg = NULL,
                                 type = c("genome", "mds"),
                                 verbose = TRUE) {
  type <- match.arg(type)

  # Check if map type is 'mds', if so, stop the process as it's not supported yet
  if(type == "mds")
    stop("Map integration is only available for genome-ordered maps. Integration for mds-ordered maps will be available soon.")

  # Parse linkage group and type for each map
  z <- lapply(x, parse_lg_and_type, lg, type)
  lg <- z[[1]]$lg

  # Ensure all elements in 'x' are of class 'mappoly2.sequence'
  assert_that(all(sapply(x, function(x) is.mappoly2.sequence(x))),
              msg = "All elements in 'x' must be of class 'mappoly2.sequence'")

  # Construct names for each biparental population
  names(x) <- sapply(x, function(x) paste0(x$data$name.p1, "_x_", x$data$name.p2))

  # Create a matrix indicating the presence of genome information for each map in 'x'
  map.mat <- lapply(x, function(x) sapply(x$maps[lg], function(x) !is.null(x[[type]][["p1p2"]])))
  temp <- Reduce(rbind, map.mat)
  rownames(temp) <- names(map.mat)
  map.mat <- t(temp)

  # Identify maps that are missing across all elements in 'x'
  y <- apply(map.mat, 1, function(x) all(!x))

  # Stop execution if there are maps missing in all populations
  if(any(y))
    stop("At least one population should have a map for group(s): ", paste(names(y)[y], collapse = " "))

  lg.names <- names(y)[!y]

  # Create a transposed matrix of parent names for each element in x
  parents.mat <- t(sapply(x, function(x) c(x$data$name.p1, x$data$name.p2)))

  # Gathering genome position for all markers
  chromo <- unique(unlist(lapply(x, function(x) unique(x$data$chrom))))
  geno.pos <- vector("list", length(chromo))
  names(geno.pos) <- chromo
  for(i in chromo){
    for(j in 1:length(x)){
      geno.pos[[i]] <- unique(rbind(geno.pos[[i]],
                                    data.frame(mrk = names(which(x[[j]]$data$chrom == i)),
                                               chr = i,
                                               pos = x[[j]]$data$genome.pos[names(which(x[[j]]$data$chrom == i))],
                                               row.names = NULL)))
    }
  }
  all.geno.pos <- NULL
  for(i in chromo)
    all.geno.pos <- rbind(all.geno.pos, geno.pos[[i]][order(geno.pos[[i]]$pos),])
  n <- all.geno.pos$mrk
  all.geno.pos <- all.geno.pos$pos
  names(all.geno.pos) <- n

  u <- vector("list", length(x))
  names(u) <- names(x)
  for(j in 1:length(x)){#population
    a1 <- reshape2::melt(x[[j]]$data$chrom)
    ## FIXME melt in a better way
    a2 <- reshape2::melt(lapply(x[[j]]$maps[lg], function(o) rownames(o$genome$p1p2$hmm.phase[[1]]$p1)))
    if(ncol(a2)==3)
      a <- cbind(a1[a2$value,,drop = FALSE], a2$Var2)
    else
      a <- cbind(a1[a2$value,,drop = FALSE], a2$L1)
    a <- unique(a)
    b <- a[,1]
    names(b) <- a[,2]
    temp <- split(names(x[[j]]$data$chrom), as.factor(x[[j]]$data$chrom))
    u[[j]] <- temp[b[names(x[[j]]$maps)]]
    u[[j]]  <- u[[j]] [!sapply(u[[j]] , is.null)]
  }
  v <- vector("list", length(lg))
  names(v) <- lg.names
  for(i in lg.names){ #lg
    for(j in 1:length(u)){
      v[[i]] <- unique(c(v[[i]], u[[j]][[b[i]]]))
    }
    v[[i]] <- sort(all.geno.pos[v[[i]]])
  }

  # Gathering parent's phases and ploidy levels
  pt <- as.vector(t(parents.mat))
  w <- table(pt)[unique(pt)]
  hom.res <- phases <- vector("list", length(w))
  names(hom.res) <- names(phases) <- names(w)

  pl.temp <- NULL
  for(i in 1:length(x)) {
    pl.temp <- rbind(pl.temp, data.frame(parent = c(x[[i]]$data$name.p1, x[[i]]$data$name.p2),
                                         ploidy = c(x[[i]]$data$ploidy.p1, x[[i]]$data$ploidy.p2)))
  }

  # Validate consistent ploidy levels across populations
  parents <- unique(pl.temp$parent)
  pl <- numeric(length(parents))
  names(pl) <- parents
  for(i in parents) {
    if(length(unique(pl.temp$ploidy[pl.temp$parent == i])) != 1)
      stop("Parent ", i, " has different ploidy levels across populations")
    pl[i] <- unique(pl.temp$ploidy[pl.temp$parent == i])
  }

  # Gathering pedigree information
  pedigree <- NULL
  for(i in 1:length(x)){
    id <- match(c(x[[i]]$data$name.p1, x[[i]]$data$name.p2), parents)
    temp <- data.frame(ind = x[[i]]$data$screened.data$ind.names,
                       Par1 = id[1],
                       Par2 = id[2],
                       pl1 = as.integer(pl[id][1]),
                       pl2 = as.integer(pl[id][2]),
                       pop = i)
    temp <- tibble::column_to_rownames(temp, var = "ind")
    pedigree <- rbind(pedigree, temp)
  }

  # Initialize results container
  homolog.correspondence <- results <- vector("list", length(lg))
  names(homolog.correspondence) <- names(results) <- lg.names

  # Process each linkage group
  for(j in lg.names){
    # Process each unique parent
    for(i in names(phases)) {
      par.ord <- which(parents.mat == i, arr.ind = TRUE)
      colnames(par.ord) <- c("pop", "parent")
      par.ord <- par.ord[order(par.ord[,1]), ,drop = FALSE]
      hom.res[[i]] <- match_homologs(x, par.ord, lg = j, pl = pl[i], par.name = i)
      ph.temp <- matrix(NA, nrow = length(v[[j]]), ncol = pl[i],
                        dimnames = list(names(v[[j]]), colnames(hom.res[[i]]$ph)))
      ph.temp[rownames(hom.res[[i]]$ph),] <- hom.res[[i]]$ph
      ### Checking dose consistency across populations
      bla <- apply(par.ord, 1, function(y) x[[y[1]]]$data[c("dosage.p1", "dosage.p2")][[y[2]]])
      ble <- unique(unlist(lapply(bla, names)))
      bli <- matrix(NA, nrow = length(ble), ncol = length(bla), dimnames = list(ble, names(bla)))
      for(k in 1:ncol(bli)){
        bli[names(bla[[k]]),k] <- bla[[k]]
      }
      untrusty.mrk <- names(which(apply(bli, 1, function(x) length(unique(x)) != 1)))
      ph.temp[intersect(rownames(ph.temp), untrusty.mrk),][]<-NA
      phases[[i]] <- ph.temp
    }
    id <- apply(sapply(phases, function(x) apply(x, 1, sum)), 1, function(x) !all(is.na(x)))
    for(i in 1:length(phases)){
      phases[[i]] <- phases[[i]][id,]
    }
    # Construct the genetic matrix
    G <- matrix(NA, nrow = nrow(phases[[1]]), ncol = nrow(pedigree),
                dimnames = list(rownames(phases[[1]]), rownames(pedigree)))

    for(i in 1:nrow(parents.mat)) {
      z <- x[[i]]$data$geno.dose
      row_intersect <- intersect(rownames(z), rownames(phases[[1]]))
      col_intersect <- intersect(colnames(z), rownames(pedigree))
      G[row_intersect, col_intersect] <- z[row_intersect, col_intersect]
    }

    # Store results for the current linkage group
    results[[j]] <- list(PH = phases, G = G, pedigree = pedigree)

    homolog.correspondence[[j]] <- hom.res
  }

  # Return results with a specific class
  structure(list(phases = results,
                 pedigree = pedigree,
                 homolog.correspondence = homolog.correspondence,
                 ploidy = pl,
                 individual.maps = x),
            class = "mappoly2.prepared.integrated.data")
}


construct_dose_matrix <- function(phases, pedigree, parents.mat, x) {
  G <- matrix(NA, nrow = nrow(phases[[1]]), ncol = nrow(pedigree),
              dimnames = list(rownames(phases[[1]]), rownames(pedigree)))

  for(i in 1:nrow(parents.mat)) {
    z <- x[[i]]$data$geno.dose
    row_intersect <- intersect(rownames(z), rownames(phases[[1]]))
    col_intersect <- intersect(colnames(z), rownames(pedigree))
    G[row_intersect, col_intersect] <- z[row_intersect, col_intersect]
  }
  return(G)
}

#' @export
print.mappoly2.prepared.integrated.data  <- function(x,...){
  y <- table(x$pedigree$Par1, x$pedigree$Par2)
  parent.names <- names(x$phases[[1]]$PH)
  dimnames(y) <- list(parent.names[as.numeric(rownames(y))], parent.names[as.numeric(colnames(y))])
  fds <- length(parent.names)
  n.mrk <- sapply(x$phases, function(x) nrow(x$PH[[1]]))
  {
    cat("Founder names:           ",  paste(parent.names, collapse = " / "), "\n")
    cat("Ploidy of founders:      ", x$ploidy, "\n")
    cat("No. individuals:         ", sum(y), "\n")
    cat("No. markers              ", sum(n.mrk), "\n\n")
    cat("Number of individuals per crosses:\n")
    print_matrix(y, spaces = 0, equal.space = FALSE)
  }
}

#' @export
plot.mappoly2.prepared.integrated.data  <- function(x, lg = 1, ...){
  pl <- x$ploidy
  assert_that(is.numeric(lg) & lg <= length(x$phases))
  w <- x$homolog.correspondence[[lg]]
  hc <- lapply(w, function(x) x$hc)
  hc <-  hc[!sapply(hc, function(x) all(is.na(x)))]
  a <- optimal_layout(length(hc))
  flag <- FALSE
  if(identical(a, c(1,1))){
    a <- c(1,3)
    flag <- TRUE
  }
  op <- par(mfrow = a, pty = "s")
  on.exit(par(op))
  for(i in 1:length(hc)){
    if(flag)
      plot(0,0, type = "n", axes = FALSE, xlab = "", ylab = "")
    d <- as.dendrogram(hc[[i]])
    d <- d %>%
      dendextend::color_branches(k = pl[i], col = mp_pal(pl[i])) %>%
      dendextend::color_labels(k = pl[i], col = mp_pal(pl[i])) %>%
      dendextend::set("branches_lwd", 4)
    plot(d, main = paste(names(pl)[i], "-",names(x$homolog.correspondence)[lg]), axes = FALSE)
  }
  # Add the overall title
  par(mfrow = c(1,1))
  l <- -5
  if(length(hc) > 3) l <- -1
  mtext("Correspondence among homologs across populations",
        side = 3, line = l, outer = TRUE)
}

#' Estimate Consensus Genetic Map
#'
#' This function estimates a consensus genetic map among different biparental
#' populations.
#'
#' @param x A list of 'mappoly2.prepared.integrated.data' objects.
#' @param err Error rate to be used in the HMM map estimation (default is 0.0).
#' @param ncpus Number of CPUs to use for parallel processing (default is 1).
#'        Automatically adjusted to not exceed the number of available cores.
#' @param verbose Logical; if TRUE, enables the printing of progress messages (default is TRUE).
#' @param detailed_verbose Logical; if TRUE, enables the printing of detailed progress messages (default is FALSE).
#' @param tol Tolerance level for the estimation convergence (default is 10e-4).
#' @param ret_H0 Logical; if TRUE, returns null hypothesis estimates (default is FALSE).
#'
#' @return A list of results from the consensus map estimation.
#' @importFrom parallel makeCluster stopCluster detectCores mclapply clusterExport
#' @importFrom assertthat assert_that
#' @export
estimate_consensus_map <- function(x,
                                   err = 0.0,
                                   ncpus = 1,
                                   verbose = TRUE,
                                   detailed_verbose = FALSE,
                                   tol = 10e-4,
                                   ret_H0 = FALSE) {
  # Validate that 'x' is a 'mappoly2.prepared.integrated.data' object
  assert_that(inherits(x, "mappoly2.prepared.integrated.data"))

  # Detect the operating system
  os_type <- Sys.info()["sysname"]

  # Adjust the number of cores to not exceed available cores
  ncpus <- min(ncpus, detectCores())

  # Conditional execution: Parallel or Serial
  if (ncpus > 1) {
    # For parallel execution
    if (os_type == "Windows") {
      # Windows OS: Use parLapply
      cl <- makeCluster(ncpus)
      on.exit(stopCluster(cl))
      clusterExport(cl, varlist = c("err", "verbose", "detailed_verbose", "tol", "ret_H0"), envir = environment())

      # Parallel execution using parLapply
      w <- parLapply(cl, x$phases, function(xi) {
        est_hmm_map_biallelic(PH = xi$PH,
                              G = xi$G,
                              pedigree = as.matrix(xi$pedigree),
                              rf = rep(0.01, nrow(xi$G) - 1),
                              err = err,
                              verbose = verbose,
                              detailed_verbose = detailed_verbose,
                              tol = tol,
                              ret_H0 = ret_H0)
      })
    } else {
      # Non-Windows OS: Use mclapply
      w <- mclapply(x$phases, function(xi) {
        est_hmm_map_biallelic(PH = xi$PH,
                              G = xi$G,
                              pedigree = as.matrix(xi$pedigree),
                              rf = rep(0.01, nrow(xi$G) - 1),
                              err = err,
                              verbose = verbose,
                              detailed_verbose = detailed_verbose,
                              tol = tol,
                              ret_H0 = ret_H0)
      }, mc.cores = ncpus)
    }
  } else {
    # For serial execution
    w <- lapply(x$phases, function(xi) {
      est_hmm_map_biallelic(PH = xi$PH,
                            G = xi$G,
                            pedigree = as.matrix(xi$pedigree),
                            rf = rep(0.01, nrow(xi$G) - 1),
                            err = err,
                            verbose = verbose,
                            detailed_verbose = detailed_verbose,
                            tol = tol,
                            ret_H0 = ret_H0)
    })
  }
  result <- vector("list", length(x$phases))
  names(result) <- names(x$phases)
  for(i in names(result)){
    result[[i]] <- list(ph = x$phases[[i]],
                        loglike = w[[i]][[1]],
                        rf = w[[i]][[2]],
                        error = err,
                        haploprob = NULL)
  }
  # Return results with specific class
  structure(list(consensus.map = result,
                 individual.maps = x$individual.maps),
            class = "mappoly2.consensus.map")
}

#' @export
print.mappoly2.consensus.map  <- function(x,...){
  y <- table(x$consensus.map[[1]]$ph$pedigree$Par1,
             x$consensus.map[[1]]$ph$pedigree$Par2)
  parent.names <- names(x$consensus.map[[1]]$ph$PH)
  dimnames(y) <- list(parent.names[as.numeric(rownames(y))], parent.names[as.numeric(colnames(y))])
  fds <- length(parent.names)
  n.mrk <- sapply(x$consensus.map, function(x) length(x$rf)+1)
  a1<-sapply(x$consensus.map, function(x) sapply(x$ph$PH, function(x) ncol(x)))
  a2<-sapply(x$consensus.map, function(x) is.null(x$haploprob))
  msg("Consensus Map:", col = "blue")
  cat("Ploidy of founders:             ", a1[,1], "\n")
  cat("Total No. individuals:          ", sum(y), "\n")
  cat("Total No. markers               ", sum(n.mrk), "\n")
  cat("Haplotype probability computed: ", ifelse(all(a2), "No", "Yes"), "\n\n")
  cat("Number of individuals per crosse:\n")
  print_matrix(y, spaces = 0, equal.space = FALSE)
  map.len <- sapply(x$consensus.map, function(x) round(sum(imf_h(x$rf)),2))
  cat("\nConsensus Map:\n---------------\n")
  R <- data.frame('LG' = names(n.mrk),
                  'Map_length_(cM)' = map.len,
                  'Markers/cM' = round(n.mrk/map.len, 3),
                  'Total mrk' = n.mrk,
                  'Max_gap' = sapply(x$consensus.map, function(x) round(max(imf_h(x$rf)),2)))
  print_matrix(R, row.names = FALSE, spaces = 0, equal.space = FALSE)
  msg("", col = "blue")
}


#' Calculate Haplotype Probabilities for Consensus Maps
#'
#' This function calculates haplotype probabilities for each linkage group in a consensus
#' map dataset of interconnected F1 populations. It utilizes parallel processing for
#' efficient computation when multiple CPU cores are available.
#'
#' @param x An object of class 'mappoly2.consensus.map', representing consensus map data
#'          for interconnected F1 populations.
#' @param ncpus The number of CPU cores to use for parallel processing, defaults to 1.
#'          If more than 1 core is specified, parallel processing is used.
#'
#' @return The input object `x` with haplotype probabilities calculated for each consensus map.
#'
#' @details The function processes the consensus map data to calculate haplotype probabilities
#'          using a biallelic model. It leverages parallel processing capabilities when
#'          multiple cores are available, significantly improving efficiency on large datasets.
#' @export
calc_consensus_haplo <- function(x,
                                 ncpus = 1) {
  # Validate that 'x' is a 'mappoly2.prepared.integrated.data' object
  assert_that(inherits(x, "mappoly2.consensus.map"))

  # Detect the operating system
  os_type <- Sys.info()["sysname"]

  # Adjust the number of cores to not exceed available cores
  ncpus <- min(ncpus, detectCores())

  # Conditional execution: Parallel or Serial
  if (ncpus > 1) {
    # For parallel execution
    if (os_type == "Windows") {
      # Windows OS: Use parLapply
      cl <- makeCluster(ncpus)
      on.exit(stopCluster(cl))
      clusterExport(cl, varlist = c("pedigree"), envir = environment())

      # Parallel execution using parLapply
      w <- parLapply(cl, x$consensus.map, function(xi) {
        u <- calc_haploprob_biallelic(PH = xi$ph$PH,
                                      G = xi$ph$G,
                                      pedigree = as.matrix(xi$ph$pedigree),
                                      rf = xi$rf,
                                      err = xi$error)
        par_col <- unlist(apply(as.matrix(xi$ph$pedigree), 1, function(x) c(rep(x[1], x[3]),
                                                                            rep(x[2], x[4])),
                                simplify = FALSE))
        u <- cbind(par_col,u)
        colnames(u)[1:3] <- c("parent", "ind", "homolog")
        u
      })
    } else {
      # Non-Windows OS: Use mclapply
      w <- mclapply(x$consensus.map, function(xi) {
        u <- calc_haploprob_biallelic(PH = xi$ph$PH,
                                      G = xi$ph$G,
                                      pedigree = as.matrix(xi$ph$pedigree),
                                      rf = xi$rf,
                                      err = xi$error)
        par_col <- unlist(apply(as.matrix(xi$ph$pedigree), 1, function(x) c(rep(x[1], x[3]),
                                                                            rep(x[2], x[4])),
                                simplify = FALSE))
        u <- cbind(par_col,u)
        colnames(u)[1:3] <- c("parent", "ind", "homolog")
        u
      }, mc.cores = ncpus)
    }
  } else {
    # For serial execution
    w <- lapply(x$consensus.map, function(xi) {
      u <- calc_haploprob_biallelic(PH = xi$ph$PH,
                                    G = xi$ph$G,
                                    pedigree = as.matrix(xi$ph$pedigree),
                                    rf = xi$rf,
                                    err = xi$error)
      par_col <- unlist(apply(as.matrix(xi$ph$pedigree), 1, function(x) c(rep(x[1], x[3]),
                                                                          rep(x[2], x[4])),
                              simplify = FALSE))
      u <- cbind(par_col,u)
      colnames(u)[1:3] <- c("parent", "ind", "homolog")
      u
    })
  }
  for(i in 1:length(x$consensus.map))
  {
    x$consensus.map[[i]]$haploprob <- w[[i]]
  }
  return(x)
}

