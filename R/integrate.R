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
#' @noRd
#' @keywords internal
match_homologs <- function(x,
                           par.ord,
                           lg,
                           pl,
                           par.name){
  # Number of full-sibs with the analyzed parent
  n <- nrow(par.ord)

  # Identify shared markers across all maps for a given linkage group and parent
  idn <- Reduce(intersect, lapply(x, function(x) rownames(x$maps[[lg]]$genome$p1p2$hmm.phase[[1]]$p1)))

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
    ph <- x[[par.ord[i,1]]]$maps[[lg]]$genome$p1p2$hmm.phase[[1]][[par.ord[i,2]]]
    colnames(ph) <- paste0("H", 1:pl, "_pop_", par.ord[i,1], "_par_",  par.ord[i,2])
    ph.list[[i]] <- ph
    ph.mat <- rbind(ph.mat, t(ph[idn,]))
    ph.names <- c(ph.names, paste0("pop_", par.ord[i,1], "_par_",  par.ord[i,2]))
  }
  names(ph.list) <- ph.names

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
  }
  else {
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
  list(ph = ph.out, equivalence = x, hc = hc, shared.mrks = idn, genome.pos = pos)
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

  # Parse linkage group and type for each map
  z <- lapply(x, parse_lg_and_type, lg, type)
  lg <- z[[1]]$lg

  # Check if map type is 'mds', if so, stop the process as it's not supported yet
  if(type == "mds")
    stop("Map integration is only available for genome-ordered maps. Integration for mds-ordered maps will be available soon.")

  # Ensure all elements in 'x' are of class 'mappoly2.sequence'
  assert_that(all(sapply(x, function(x) is.mappoly2.sequence(x))),
              msg = "All elements in 'x' must be of class 'mappoly2.sequence'")

  # Construct names for each biparental population
  names(x) <- sapply(x, function(x) paste0(x$data$name.p1, "x", x$data$name.p2))

  # Create a matrix indicating the presence of genome information for each map in 'x'
  map.mat <- sapply(x, function(x) sapply(x$maps, function(x) !is.null(x[[type]][["p1p2"]])), simplify = "array")

  # Identify maps that are missing across all elements in 'x'
  y <- apply(map.mat, 1, function(x) all(!x))

  # Stop execution if there are maps missing in all populations
  if(any(y))
    stop("At least one population should have a map for group(s): ", paste(names(y)[y], collapse = " "))

  # Create a transposed matrix of parent names for each element in x
  parents.mat <- t(sapply(x, function(x) c(x$data$name.p1, x$data$name.p2)))

  # Gathering parent's phases and ploidy levels
  w <- table(as.vector(parents.mat))
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
  pedigree <- create_pedigree(parents.mat, x, pl)

  # Initialize results container
  results <- vector("list", length(lg))

  # Process each linkage group
  for(j in lg){
    # Process each unique parent
    for(i in names(phases)) {
      par.ord <- which(parents.mat == i, arr.ind = TRUE)
      colnames(par.ord) <- c("pop", "parent")
      par.ord <- par.ord[order(par.ord[,1]),]
      hom.res[[i]] <- match_homologs(x, par.ord, lg = j, pl = pl[i], par.name = i)
      phases[[i]] <- hom.res[[i]]$ph
    }

    # Construct the genetic matrix
    G <- construct_dose_matrix(phases, pedigree, parents.mat, x)

    # Store results for the current linkage group
    results[[j]] <- list(PH = phases, G = G, pedigree = pedigree, homolog.correspondence = hom.res, ploidy = pl)
  }

  # Return results with a specific class
  structure(results, class = "mappoly2.prepared.integrated.data")
}

create_pedigree <- function(parents.mat, x, pl) {
  pedigree <- NULL
  char_vector <- as.vector(parents.mat)
  unique_strings <- unique(char_vector)
  string_to_int <- setNames(seq_along(unique_strings), unique_strings)
  int_vector <- string_to_int[char_vector]
  par.idx <- matrix(int_vector, nrow = nrow(parents.mat))

  for(i in 1:length(x)) {
    all.ind <- x[[i]]$data$screened.data$ind.names
    temp_pedigree <- data.frame(Ind = all.ind,
                                Par1 = par.idx[i, 1],
                                Par2 = par.idx[i, 2],
                                pl1 = as.integer(pl[par.idx[i, 1]]),
                                pl2 = as.integer(pl[par.idx[i, 2]]),
                                pop = i)
    pedigree <- rbind(pedigree, temp_pedigree)
  }

  pedigree <- tibble::column_to_rownames(pedigree, var = "Ind")
  return(pedigree)
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
plot.mappoly2.prepared.integrated.data  <- function(x, lg = 1, ...){
  pl <- x[[1]]$ploidy
  assert_that(is.numeric(lg), lg <= length(x[[1]]))
  a <- mappoly2:::optimal_layout(length(pl))
  op <- par(mfrow = a, pty = "s")
  on.exit(par(op))
  k <- lg
  w <- x[[k]]
  hc <- lapply(w$homolog.correspondence, function(x) x$hc)
  hc <-  hc[!sapply(hc, function(x) all(is.na(x)))]
  for(i in 1:length(hc)){
    d <- as.dendrogram(hc[[i]])
    d <- d %>%
      dendextend::color_branches(k = pl[i], col = drsimonj_colors(pl[i])) %>%
      dendextend::color_labels(k = pl[i], col = drsimonj_colors(pl[i])) %>%
      set("branches_lwd", 4)
      plot(d, main = names(pl)[i], axes = FALSE)
  }
  # Add the overall title
  par(mfrow = c(1,1))
  mtext("Correspondence among homologs across populations",
        side = 3, line = -5, outer = TRUE)
}

