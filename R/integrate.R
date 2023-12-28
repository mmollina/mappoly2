#' Match Homologs Across Biparental Genetic Maps
#'
#' This internal function matches shared parents across biparental genetic maps and rearranges
#' the homologs to ensure consistency in joint map estimation. If inconsistencies are found,
#' such as different dosages or phasing, the function removes the marker. It operates on a single
#' chromosome.
#'
#' @param map.list A list of biparental genetic maps.
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
match_homologs <- function(map.list,
                           par.ord,
                           lg,
                           pl,
                           par.name){
  # Number of full-sibs with the analyzed parent
  n <- nrow(par.ord)

  # Identify shared markers across all maps for a given linkage group and parent
  idn <- Reduce(intersect, lapply(map.list, function(x) rownames(x$maps[[lg]]$genome$p1p2$hmm.phase[[1]]$p1)))

  # Collect all markers and their positions from each map
  pos.all <- lapply(map.list, function(x) data.frame(x$data$mrk.names, x$data$genome.pos))
  pos.all <- Reduce(rbind, pos.all)
  pos.all <- pos.all[!duplicated(pos.all[,1]),]
  dimnames(pos.all) <- list(pos.all[,1], c("mrk", "pos"))
  pos.all <- pos.all[unique(unlist(lapply(map.list, function(x) rownames(x$maps[[lg]]$genome$p1p2$hmm.phase[[1]]$p1)), use.names = FALSE)),,drop = FALSE]

  # Order markers according to genome positions
  pos <- pos.all[order(pos.all$pos),-1, drop = FALSE]

  # Collect phases for each map and parent
  ph.list <- vector("list", n)
  ph.mat <- NULL
  for(i in 1:n){
    ph <- map.list[[par.ord[i,1]]]$maps[[lg]]$genome$p1p2$hmm.phase[[1]][[par.ord[i,2]]]
    colnames(ph) <- paste0(letters[1:pl], i)
    ph.list[[i]] <- ph
    ph.mat <- rbind(ph.mat, t(ph[idn,]))
  }

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
    x <- split(names(homologs), as.factor(homologs))
    names(x) <- paste0("h", 1:pl)
    x <- x %>%
      melt %>%
      mutate(homolog = substr(value, 1,1)) %>%
      mutate(pop = substr(value,2,10)) %>%
      arrange(pop, L1) %>%
      acast(pop ~ L1, value.var = "value")

    # Reorganize homologs for consistency
    for(i in 1:length(ph.list))
      ph.list[[i]] <- ph.list[[i]][,x[i,]]

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
    for(i in 1:length(map.list)){
      idtemp <- setdiff(rownames(map.list[[i]]$maps[[lg]]$genome$p1p2$hmm.phase[[1]]$p1), rownames(ph.out))
      ph.out <- rbind(ph.out, matrix(NA, length(idtemp), pl, dimnames = list(idtemp, NULL)))
    }
  }

  # Rearrange the output according to genome positions
  ph.out <- ph.out[rownames(pos),]
  colnames(pos) <- "geno.pos"

  # Return the final list of outputs
  list(ph = ph.out, equivalence = x, hc = hc, shared.mrks = idn, genome.pos = pos)
}

#' @export
#' @importFrom ggplot2 ggplot aes geom_point xlab ylab theme
plot_multi_map <- function(x){
  assert_that(all(sapply(x, function(x) is.mappoly2.sequence(x))),
              msg = "all elements in 'x' must be of class 'mappoly2.sequence'")
  names(x) <- sapply(x, function(x) paste0(x$data$name.p1 , "x", x$data$name.p2))
  map.mat <- sapply(x, function(x)  sapply(x$maps, function(x) !is.null(x[["genome"]][["p1p2"]])), simplify = "array")
  map.mat <- matrix(map.mat, nrow = length(x[[1]]$maps), dimnames = list(names(x[[1]]$maps), names(x)))
  y <- apply(map.mat, 1, function(x) all(!x))
  if(any(y))
    stop("at least one population should have a map for group(s): ", paste(names(y)[y], collapse = " "))
  maps <- NULL

  for(i in 1:length(x)){
    w <- x[[i]]
    for(j in rownames(map.mat)[map.mat[,i]]){
      maps <- rbind(maps, data.frame(F1 = names(x)[i],
                                     LG = j,
                                     mrk.names = rownames(w$maps[[j]]$genome$p1p2$hmm.phase[[1]]$p1),
                                     pos = cumsum(c(0, imf_h(w$maps[[j]]$genome$p1p2$hmm.phase[[1]]$rf)))))
    }
  }
  # Count the occurrences of each marker name in 'mrk.names'
  marker_counts <- maps %>%
    count(mrk.names) %>%
    mutate(category = as.character(n))

  # Join the counts back to the original data frame
  maps_with_counts <- maps %>%
    left_join(marker_counts, by = "mrk.names")

  # Generate a color set for each unique count category using viridis
  unique_counts <- sort(unique(marker_counts$category))
  colors <- rev(viridis::mako(length(unique_counts)))
  names(colors) <- unique_counts

  # Update the plotting code
  lo <- optimal_layout(nrow(map.mat))
  ggplot(maps_with_counts, aes(x = pos, y = F1, group = as.factor(F1), color = category)) +
    geom_point(shape = 108, size = 5, show.legend = TRUE) +
    facet_wrap(vars(LG), nrow = lo[1], ncol = lo[2]) +
    xlab("Position (cM)") +
    ylab("Biparental Maps") +
    scale_color_manual(values = colors, name = "Marker\nFrequency\nAcross\nPopulations") +  # Set the legend title
    theme(legend.title = element_text(face = "bold")) +  # Optionally make the legend title bold
    guides(color = guide_legend(override.aes = list(shape = 15))) + # Square shape for legend keys
    theme_dark()
}

#' @export
prepare_to_integrate <- function(map.list,
                                 lg) {

  # Ensure that all elements in 'x' are of class 'mappoly2.sequence'
  assert_that(all(sapply(map.list, function(x) is.mappoly2.sequence(x))),
              msg = "all elements in 'x' must be of class 'mappoly2.sequence'")

  # Rename each element in 'x' based on the concatenation of parent names
  names(map.list) <- sapply(map.list, function(x) paste0(x$data$name.p1 , "x", x$data$name.p2))

  # Create a matrix (map.mat) indicating the presence of genome information for each map in 'x'
  map.mat <- sapply(map.list, function(x) sapply(x$maps, function(x) !is.null(x[["genome"]][["p1p2"]])), simplify = "array")

  # Reshape map.mat into a 2D matrix with rows as maps and columns as elements in 'x'
  map.mat <- matrix(map.mat, nrow = length(map.list[[1]]$maps), dimnames = list(names(map.list[[1]]$maps), names(x)))

  # Identify maps that are missing across all elements in 'x'
  y <- apply(map.mat, 1, function(x) all(!x))

  # Stop execution if there are maps missing in all populations
  if(any(y))
    stop("at least one population should have a map for group(s): ", paste(names(y)[y], collapse = " "))

  # Create a transposed matrix of parent names for each element in map.list
  parents.mat <- t(sapply(map.list, function(x) c(x$data$name.p1, x$data$name.p2)))

  # Set names for each map in the list based on parent names
  names(map.list) <- sapply(map.list, function(x) paste0(x$data$name.p1, "x", x$data$name.p2))

  # Gathering parent's phases
  w <- table(as.vector(parents.mat))
  hom.res <- phases <- vector("list", length(w))
  names(hom.res) <- names(phases) <- names(w)

  # Gathering ploidy levels
  pl.temp <- NULL
  for(i in 1:length(map.list)) {
    pl.temp <- rbind(pl.temp, data.frame(parent = c(map.list[[i]]$data$name.p1,map.list[[i]]$data$name.p2),
                                         ploidy = c(map.list[[i]]$data$ploidy.p1, map.list[[i]]$data$ploidy.p2)))
  }
  parents <- unique(pl.temp$parent)
  pl <- numeric(length(parents))
  names(pl) <- parents
  for(i in parents) {
    if(length(unique(pl.temp$ploidy[pl.temp$parent == i])) == 1)
      pl[i] <- unique(pl.temp$ploidy[pl.temp$parent == i])
    else
      stop("parent ", i, " has different ploidy levels across populations")
  }

  # Process each unique parent
  for(i in names(phases)) {
    par.ord <- which(parents.mat == i, arr.ind = TRUE)
    colnames(par.ord) <- c("pop", "parent")
    hom.res[[i]] <- match_homologs(map.list, par.ord, lg, pl = pl[i], par.name = i)
    phases[[i]] <- hom.res[[i]]$ph
  }

  # Gathering pedigree information
  pedigree <- NULL
  char_vector <- as.vector(parents.mat)
  unique_strings <- unique(char_vector)
  string_to_int <- setNames(seq_along(unique_strings), unique_strings)
  int_vector <- string_to_int[char_vector]
  par.idx <- matrix(int_vector, nrow = nrow(parents.mat))
  all.ind <- lapply(map.list, function(x) x$data$screened.data$ind.names)
  all.ind.id <- NULL
  for(i in 1:length(all.ind)) {
    all.ind.id <- c(all.ind.id, all.ind[[i]])
    pedigree <- rbind(pedigree, data.frame(Ind = all.ind[[i]],
                                           Par1 = par.idx[i,1],
                                           Par2 = par.idx[i,2],
                                           pl1 = as.integer(pl[par.idx[i,1]]),
                                           pl2 = as.integer(pl[par.idx[i,2]]),
                                           pop = i))
  }
  pedigree <- column_to_rownames(pedigree, "Ind")

  # Construct the genetic matrix
  G <- matrix(NA, nrow(phases[[1]]), nrow(pedigree), dimnames = list(rownames(phases[[1]]), rownames(pedigree)))
  for(i in 1:nrow(parents.mat)) {
    z <- map.list[[i]]$data$geno.dose
    G[intersect(rownames(z), rownames(phases[[1]])), intersect(colnames(z), rownames(pedigree))] <-
      z[intersect(rownames(z), rownames(phases[[1]])), intersect(colnames(z), rownames(pedigree))]
  }

  # Return a list of results
  list(PH = phases, G = G, pedigree = pedigree, homolog.correspondence = hom.res)
}

#' @export
plot_prepared_data <- function(w){
  hc <- lapply(w$homolog.correspondence, function(x) x$hc)
  hc <-  hc[!sapply(hc, function(x) all(is.na(x)))]
  #u1 <- unique(w$pedigree[,c(1,3)])
  #v1 <- u1[,2]
  #names(v1) <- u1[,1]
  #u2 <- unique(w$pedigree[,c(2,4)])
  #v2 <- u2[,2]
  #names(v2) <- u2[,1]
  #v <- c(v1, v2)
  #names(v) <- names(hc)
  a <- mappoly2:::optimal_layout(length(hc))
  par(mfrow = a)
  for(i in 1:length(hc)){
    d <- as.dendrogram(hc[[i]])
    #d <- d %>%
    #  dendextend::color_branches(k = v[names(hc)[1]], col = mappoly::mp_pallet2(v[names(hc)[1]])) %>%
    #  dendextend::color_labels(k = v[names(hc)[1]], col = mappoly::mp_pallet2(v[names(hc)[1]]))
    plot(d, main = names(hc)[i])
    #dendextend::rect.dendrogram(d, k = v[names(hc)[1]], lwd = 3,
    #                            border = mappoly::mp_pallet2(v[names(hc)[1]]))
  }
}

