#' Given a screened data set, returns the indices of the
#' screened markers in the corresponding raw dataset
#' @keywords internal
get_screened_mrk_indices <- function(x){
  match(x$screened.data$mrk.names, x$mrk.names)
}

#' Given a dataset containing chromosome information,
#' and a chromosome vector return a vector of screened
#' markers corresponding to the chromosomes provided
#' @keywords internal
get_mrk_indices_from_chrom <- function(x, chrom){
  assert_that(has.chromosome.info(x))
  a <- unique(x$chrom)
  na.name <- a[is.na(embedded_to_numeric(a))]
  ch.n.arg <- embedded_to_numeric(chrom)
  if(any(is.na(ch.n.arg))){
    b <- chrom%in%na.name & is.na(ch.n.arg)
    if(all(!b))
      stop(msg = "provide a valid chromosome identifier")
  }
  ch.n.dat <- embedded_to_numeric(x$chrom)
  if(any(is.na(ch.n.arg)))
    ch.n.arg[is.na(ch.n.arg)] <- ch.n.dat[is.na(ch.n.dat)] <- 999
  ch.id <- which(ch.n.dat%in%ch.n.arg)
  if(has.mappoly2.screened(x))
    ch.id <- intersect(ch.id, get_screened_mrk_indices(x))
  ch.id
}

#' Removes all non-digit characters and then converts
#' the remaining string to an integer.
#' @param void internal function
#' @keywords internal
embedded_to_numeric <- function(x) {
  as.integer(gsub("[^0-9]", "", x))
}

#' Creates a formatted message with a horizontal rule
#' and specified text color.
#' @importFrom cli rule
#' @importFrom magrittr %>%
#' @keywords internal
msg <- function(text, line = 1, col = "-"){
  cli::rule(line = line, left = text) %>%
    text_col(col = col) %>%
    message()
}

#' Changes the color of text in the R console, with a specific
#' focus on compatibility with RStudio's dark and light themes.
#' @importFrom rstudioapi isAvailable hasFun getThemeInfo
#' @importFrom crayon white black
#' @keywords internal
text_col <- function(x, col = c("-","red", "blue")) {
  col <- match.arg(col)
  if(col == "red") return(crayon::red(x))
  else if(col == "blue") return(crayon::blue(x))
  # If RStudio not available, messages already printed in black
  if (!rstudioapi::isAvailable()) {
    return(x)
  }
  if (!rstudioapi::hasFun("getThemeInfo")) {
    return(x)
  }
  theme <- rstudioapi::getThemeInfo()
  if (isTRUE(theme$dark)) crayon::white(x) else crayon::black(x)
}

#' Sort Phases Based on Log-Likelihood
#'
#' This function sorts a list of phases based on their log-likelihood values.
#' It can return either the best phase (highest log-likelihood) or all phases in a sorted order.
#'
#' @param ph A list of phases, each phase should have a `$loglike` element.
#' @param only.best A boolean flag; if TRUE, only the best (highest log-likelihood) phase is returned.
#'                   If FALSE, all phases are returned in descending order of their log-likelihoods. Default is TRUE.
#' @return If `only.best` is TRUE, returns the phase with the highest log-likelihood.
#'         If FALSE, returns all phases sorted by log-likelihood in descending order.
#' @keywords internal
sort_phase <- function(ph, only.best = TRUE){
  ph <- ph[sapply(ph, function(x) !is.null(x$loglike))]
  if(only.best){
    ph <- ph[which.max(sapply(ph, function(x) x$loglike))]
    return(ph)
  } else {
    ph <- ph[order(sapply(ph, function(x) x$loglike), decreasing = TRUE)]
    return(ph)
  }
}

#' Find Flanking Markers
#'
#' Identifies the closest preceding and succeeding elements from `B1` for each element in `B2`,
#' based on their positions in `A`.
#'
#' @param A A character vector; positions of elements in this vector are considered.
#' @param B1 A character vector; contains reference elements.
#' @param B2 A character vector; elements for which closest flanking elements in `B1` are to be found.
#' @return A list where each element in `B2` is mapped to its closest preceding and succeeding elements from `B1`.
#' @keywords internal internal
find_flanking_markers <- function(A, B1, B2) {
  A <- as.character(A)
  B1 <- as.character(B1)
  B2 <- as.character(B2)

  # Match B1 elements to their positions in A
  B1_positions <- match(B1, A)
  names(B1_positions) <- B1

  flanking_letters <- lapply(B2, function(b) {
    b_position <- match(b, A)
    preceding <- B1_positions[B1_positions < b_position]
    succeeding <- B1_positions[B1_positions > b_position]

    # Get the closest preceding and succeeding elements from B1
    preceding <- if (length(preceding) > 0) names(preceding)[which.max(preceding)] else NA
    succeeding <- if (length(succeeding) > 0) names(succeeding)[which.min(succeeding)] else NA

    list(preceding = preceding, succeeding = succeeding)
  })

  names(flanking_letters) <- B2
  flanking_letters
}



#' Parse Linkage Group and Type
#'
#' This function parses the linkage group and type from a given mappoly2.sequence object.
#' It validates the input and extracts the specified linkage group and type.
#'
#' @param x An object of class mappoly2.sequence. This object should contain SNP data
#'   along with associated maps.
#' @param lg A character or numeric vector specifying the linkage groups. If a character
#'   vector, it should correspond to the names of the maps in the mappoly2.sequence object.
#'   If numeric, it should be the indices of the maps.
#' @param type A character vector specifying the type of SNP order.
#'   Possible values are "mds", "genome", or "custom". Default is c("mds", "genome", "custom").
#'   The function will match the provided type argument with these options.
#'
#' @return A list containing two elements: 'lg' and 'type'. The 'lg' element contains the
#'   linkage group information, either as a character or numeric vector, as provided by the user.
#'   The 'type' element contains the matched SNP order type.
#'
#' @importFrom assertthat assert_that
#' @keywords internal
parse_lg_and_type <- function(x, lg = NULL, type = c("mds", "genome", "custom")){
  assert_that(is.mappoly2.sequence(x))
  type <- match.arg(type)
  if(is.null(lg))
    lg <- seq_along(x$maps)
  if(is.character(lg)){
    assert_that(lg %in% names(x$maps), msg = "Provide a valid group set")
    lg <- match(lg, names(x$maps))
  }
  assert_that(all(lg <= length(x$maps)))
  return(list(lg = lg, type = type))
}


get_info_markers <- function(mrk.names, x, parent = c("p1p2","p1","p2")){
  parent <- match.arg(parent)
  mrk.names <- intersect(mrk.names, x$data$screened.data$mrk.names)
  if(parent == "p1")
    return(mrk.names[x$data$dosage.p2[mrk.names] == 0 | x$data$dosage.p2[mrk.names] == x$data$ploidy.p2])
  else if(parent == "p2")
    return(mrk.names[x$data$dosage.p1[mrk.names] == 0 | x$data$dosage.p1[mrk.names] == x$data$ploidy.p1])
  else if (parent == "p1p2")
    return(mrk.names)
}

detect_info_parent <- function(mrk.names, x){
  if(all(x$data$dosage.p2[mrk.names] == 0 | x$data$dosage.p2[mrk.names] == x$data$ploidy.p2))
    return("p1")
  else if(all(x$data$dosage.p1[mrk.names] == 0 | x$data$dosage.p1[mrk.names] == x$data$ploidy.p1))
    return("p2")
  else
    return("p1p2")
}

#' Extract Markers from Ordered Sequence
#'
#' Extracts markers from an ordered sequence in an object, based on the specified type.
#'
#' @param x An object containing marker maps.
#' @param lg Either a character vector naming the maps or a numeric vector indicating map indices.
#' @param type The type of sequence to consider (either "mds", "genome", or "custom").
#' @return Markers based on the specified type.
#' @keywords internal
get_markers_from_ordered_sequence <- function(x, lg, type = c("mds", "genome", "custom"),
                                              parent = c("p1p2","p1","p2")){
  y <- parse_lg_and_type(x,lg,type)
  parent <- match.arg(parent)
  if(y$type == "mds")
    w <-lapply(x$maps[y$lg], function(z, type) intersect(z[[type]]$order$locimap$locus,
                                                         z[[type]]$mkr.names), y$type)
  else if(y$type == "genome")
    w <- lapply(x$maps[y$lg], function(z, type) intersect(rownames(z[[type]]$order),
                                                          z[[type]]$mkr.names), y$type)
  w <- lapply(w, get_info_markers, x, parent)
  return(w)
}

#' Extract Markers from Mapped Sequence
#'
#' Similar to `get_markers_from_ordered_sequence` but specifically for mapped sequences.
#'
#' @param x An object containing marker maps.
#' @param lg Either a character vector naming the maps or a numeric vector indicating map indices.
#' @param type The type of sequence to consider (either "mds", "genome", or "work").
#' @return Markers based on the specified type for mapped sequences.
#' @keywords internal
get_markers_from_phased_sequence <- function(x, lg,
                                             type = c("mds", "genome", "custom"),
                                             parent = c("p1p2","p1","p2"),
                                             phase = c("rf.phase", "hmm.phase")){
  y <- parse_lg_and_type(x,lg,type)
  parent <- match.arg(parent)
  return(lapply(x$maps[y$lg], function(z, type) intersect(rownames(z[[type]][[parent]][[phase]][[1]]$p1),
                                                          z[[type]]$mkr.names), y$type))


}

#' Determine Dosage Type of Markers
#'
#' Classifies each marker based on its dosage type using the provided dosage data.
#'
#' @param x An object containing dosage data.
#' @param mrk.names Names of the markers to be analyzed.
#' @return A list categorizing each marker as simplex in parent 1, simplex in parent 2,
#'         double simplex, or multiplex based on their dosage values.
#' @keywords internal
get_dosage_type <- function(x, mrk.names){
  p1 <- abs(abs(x$data$dosage.p1 - x$data$ploidy.p1/2) - x$data$ploidy.p1/2)
  p2 <- abs(abs(x$data$dosage.p2 - x$data$ploidy.p2/2) - x$data$ploidy.p2/2)
  s.p1 <- p1  ==  1 & p2  ==  0
  s.p2 <- p1  ==  0 & p2  ==  1
  ds <- p1  ==  1 & p2  ==  1
  list(simplex.p1 = names(which(s.p1[mrk.names])),
       simplex.p2 = names(which(s.p2[mrk.names])),
       double.simplex = names(which(ds[mrk.names])),
       multiplex = names(which(!(s.p1 | s.p2 | ds)[mrk.names])))
}

#' Create a Map Skeleton for mappoly2 Data
#'
#' This function creates a skeleton structure for mapping data, specifically for mappoly2 data objects.
#' It initializes map structures for multiple linkage groups based on provided marker IDs.
#'
#' @param x A mappoly2 data object; the function checks if the input is a valid mappoly2 data structure.
#' @param mrk.id.list A list or a vector of marker IDs for which the map skeleton is to be created.
#'                    If a vector is provided, it is converted into a list.
#' @return Returns a mappoly2 sequence object with initialized map structures for each linkage group.
#'         Each linkage group will have a map for MDS, genome, and work, but with `NULL` values for `order` and `phase`.
#' @note This function is internal and usually not called directly by the user.
#' @keywords internal
.map_skeleton<- function(x, mrk.id.list){
  assert_that(is.mappoly2.data(x))
  if(!is.list(mrk.id.list))
    mrk.id.list <- list(mrk.id.list)
  n.lg <- length(mrk.id.list)
  maps <- vector("list", n.lg)
  names(maps) <- c(paste0("lg", 1:n.lg))
  for(i in 1:n.lg){
    maps[[i]] <- list(mds = list(mkr.names = mrk.id.list[[i]],
                                 order = NULL,
                                 p1 = list(rf.phase = NULL, hmm.phase = NULL),
                                 p2 = list(rf.phase = NULL, hmm.phase = NULL),
                                 p1p2 = list(rf.phase = NULL, hmm.phase = NULL)),
                      genome = list(mkr.names = mrk.id.list[[i]],
                                    order = NULL,
                                    p1 = list(rf.phase = NULL, hmm.phase = NULL),
                                    p2 = list(rf.phase = NULL, hmm.phase = NULL),
                                    p1p2 = list(rf.phase = NULL, hmm.phase = NULL)),
                      custom = list(mkr.names = mrk.id.list[[i]],
                                    order = NULL,
                                    p1 = list(rf.phase = NULL, hmm.phase = NULL),
                                    p2 = list(rf.phase = NULL, hmm.phase = NULL),
                                    p1p2 = list(rf.phase = NULL, hmm.phase = NULL)))
  }
  structure(list(maps = maps, data = x), class = "mappoly2.sequence")
}

#' Set Quality Assurance and Quality Control (QA/QC) Values
#'
#' This internal function creates a list of QA/QC values for genetic markers and individuals.
#' It constructs two data frames: one for markers and another for individuals, each containing relevant QA/QC metrics.
#'
#' @param id.mrk A vector of marker identifiers.
#' @param id.ind A vector of individual identifiers.
#' @param miss.mrk A numeric vector representing the missing data rate for each marker (default is NA).
#' @param miss.ind A numeric vector representing the missing data rate for each individual (default is a vector of NA of length equal to `id.ind`).
#' @param chisq.pval A numeric vector of chi-squared test p-values for each marker.
#' @return A list with two elements:
#'   - `markers`: A data frame with the missing data rate (`miss`), chi-squared test p-values (`chisq.pval`), and read depth (set to `NA`) for each marker.
#'   - `individuals`: A data frame with the missing data rate (`miss`) and full sibling status (`full.sib`, set to `NA`) for each individual.
#' @keywords internal
.setQAQC <- function(id.mrk, id.ind,
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

#' Retrieve Markers and Individuals Based on QA/QC Thresholds
#'
#' This internal function selects markers and individuals based on given QA/QC thresholds. It returns a list containing the names of markers and individuals that meet the specified criteria.
#'
#' @param x A list containing QA/QC values for markers and individuals, typically generated by \code{\link{.setQAQC}}.
#' @param miss.mrk.thresh A numeric value specifying the threshold for the missing data rate in markers (default is +Inf).
#' @param miss.ind.thresh A numeric value specifying the threshold for the missing data rate in individuals (default is +Inf).
#' @param chisq.pval.thresh A numeric value specifying the threshold for the chi-squared test p-value in markers (default is -Inf).
#' @param read.depth.thresh A numeric vector of length two specifying the lower and upper bounds for read depth in markers (default is c(0, +Inf)).
#' @return A list with three elements:
#'   - `thresholds`: A list of the thresholds used for selection.
#'   - `mrk.names`: The names of the markers that meet the specified criteria.
#'   - `ind.names`: The names of the individuals that meet the specified criteria.
#' @keywords internal
.get_mrk_ind_from_QAQC <- function(x,
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


#' Reverse the Order of a Genetic Map
#'
#' This function reverses the order of a genetic map in a `mappoly` object.
#' It is typically used to align the MDS order with the genomic order.
#'
#' @param x A `mappoly` object containing the genetic map(s) to be reversed.
#' @param lg An integer or vector specifying the linkage group(s) to be reversed.
#' @param type Character vector indicating the type of genetic map to be reversed.
#'             Options include "mds", "genome", or "custom". Default is c("mds", "genome", "custom").
#' @param parent A character vector specifying the parent or parents to be considered
#'               in the mapping process. Options are "p1p2" (both parents), "p1" (first parent),
#'               and "p2" (second parent). Default is c("p1p2", "p1", "p2").
#'
#' @return Returns the `mappoly` object with the order of the specified linkage group(s) reversed.
#'
#' @details The function modifies the order of markers in the specified linkage group(s) of
#'          the genetic map. It reverses the positions of markers along with associated
#'          data such as recombination fractions and haplotype probabilities, ensuring
#'          that the genetic map's orientation is consistent with genomic data.
#'
#' @export
rev_map <- function(x, lg,
                    type = c("mds", "genome", "custom"),
                    parent = c("p1p2","p1","p2")){
  y <- parse_lg_and_type(x,lg,type)
  parent <- match.arg(parent)
  for(i in y$lg){
    z<- x$maps[[i]][[y$type]][[parent]]$hmm.phase[[1]]
    z$p1 <- z$p1[nrow(z$p1):1,]
    z$p2 <- z$p2[nrow(z$p2):1,]
    z$rf <- rev(z$rf)
    if(!is.null(z$haploprob))
      z$haploprob <- cbind(z$haploprob[,1:3], z$haploprob[,ncol(z$haploprob):4])
    x$maps[[i]][[y$type]][[parent]]$hmm.phase[[1]] <- z
  }
  return(x)
}

#' Convert MAPpoly Data to CSV
#'
#' This function takes MAPpoly data objects and writes them to CSV files. Each object in the list
#' is processed individually, and a CSV file is generated for each. The user can specify the path
#' for saving the CSV files, and custom names for the parents in the dataset.
#'
#' @param x A list of MAPpoly data objects or a single MAPpoly data object. The function
#'   checks if `x` is a list, and if not, it converts it into a list.
#' @param path Optional; a string specifying the directory where the CSV files will be saved.
#'   If not provided, files are saved in the current working directory.
#' @param parent.names Optional; a matrix of names for the parent genotypes. If not provided,
#'   default names are generated in the format P1, P2, etc. The matrix should have two columns
#'   and a number of rows equal to the length of `x`.
#' @details This function processes each MAPpoly data object, extracting relevant data for
#'   CSV output. It handles missing reference and alternate sequences, and ensures that the
#'   provided parent names match the structure of the data. The function creates one CSV file
#'   per MAPpoly data object.
#' @return Invisible NULL. The function is used for its side effect of writing files.
#' @importFrom utils write.csv
#' @export
mappoly_to_csv_mappoly2 <- function(x, path = NULL, parent.names = NULL){
  if(!is.list(x) & inherits(x, "mappoly.data"))
    x <- list(x)
  if(is.null(parent.names))
    parent.names <- t(apply(matrix(1:(length(x) * 2), ncol = 2, byrow = TRUE),
                            1, function(x) paste0("P",x)))
  assert_that(is.matrix(parent.names))
  assert_that(ncol(parent.names)==2 & nrow(parent.names)==length(x))
  assert_that(all(sapply(x, function(x) inherits(x, "mappoly.data"))))
  for(i in 1:length(x)){
    F1<-x[[i]]$geno.dose
    F1[F1==x[[i]]$ploidy+1] <- NA
    r <- is.null(x[[i]]$seq.ref)
    if(is.null(r))
      r <- rep(NA, length(x[[i]]$mrk.names))
    a <- is.null(x[[i]]$seq.alt)
    if(is.null(r))
      a <- rep(NA, length(x[[i]]$mrk.names))
    w<-data.frame(snp_id = x[[i]]$mrk.names,
                  P1 = x[[i]]$dosage.p1,
                  P2 = x[[i]]$dosage.p2,
                  chrom = x[[i]]$chrom[x[[i]]$mrk.names],
                  genome_pos = x[[i]]$genome.pos[x[[i]]$mrk.names],
                  ref = r,
                  alt = a,
                  F1)
    names(w)[c(2:3)] <- parent.names[i,]
    fn <- paste0(paste0(parent.names[i, ], collapse = "x"), ".csv")
    if (is.null(path))
      file_path <- file.path(getwd(), fn)
    else
      file_path <- file.path(path, fn)
    write.csv(x, file = file_path, row.names = FALSE)
  }
}




#' MAPpoly v1 main palette
#'
#' @param void internal function to be documented
#' @keywords internal
#' @export
mp_pal <- function(n) {
  color_set <- c("#ffe119",#1
                 "#f58231",#2
                 "#e6194b",#3
                 "#800000",#4
                 "#911eb4",#5
                 "#0202a1",#6
                 "#4363d8",#7
                 "#42d4f4",#8
                 "#469990",#9
                 "#3cb44b")#10
  if (n == 2) {
      return(color_set[c(3,7)])
  } else if(n ==3){
    color_set[c(3,7,10)]
  } else if(n ==4){
    color_set[c(3,7,10,5)]
  } else if(n ==5){
    color_set[c(5,8,10,1,3)]
  } else if(n ==6){
    color_set[c(3,2,1,10,9,6)]
  } else {
    color_palette <- colorRampPalette(color_set)
    return(color_palette(n))
  }
}


detect_hmm_est_map <- function(x){
  assert_that(is.mappoly2.sequence(x))
  z <- lapply(x$maps, function(x) sapply(x, function(x) sapply(x[3:5], function(x) !is.null(x$hmm.phase))))
  v <- unlist(z)
  dim(v) <- c(3,3,length(z))
  dimnames(v) = list(c("p1", "p2", "p1p2"), c("mds", "genome", "both"), names(z))
  for(i in 1:dim(v)[3])
      v[,3,i] <- apply(v[,1:2,i],1, all)
  v
}


detect_comp_haplotype <- function(x){
  assert_that(is.mappoly2.sequence(x))
  z <- lapply(x$maps, function(x) sapply(x, function(x) sapply(x[3:5], function(x) !is.null(x$hmm.phase[[1]]$haploprob))))
  v <- unlist(z)
  dim(v) <- c(3,3,length(z))
  dimnames(v) = list(c("p1", "p2", "p1p2"), c("mds", "genome", "both"), names(z))
  for(i in 1:dim(v)[3])
    v[,3,i] <- apply(v[,1:2,i],1, all)
  v
}


get_palette <- function(all_parents, parent, n) {
  if (parent == levels(all_parents)[1])
    return(colorRampPalette(c("lightblue", "darkblue"))(n))
  else if (parent == levels(all_parents)[2])
    return(colorRampPalette(c("lightcoral", "darkred"))(n))
  else if (parent == levels(all_parents)[3])
    return(colorRampPalette(c("lightgreen", "darkgreen"))(n))
  else if (parent == levels(all_parents)[4])
    return(colorRampPalette(c("gold", "goldenrod2"))(n))
  else if (parent == levels(all_parents)[5])
    return(colorRampPalette(c("mediumpurple", "mediumpurple4"))(n))
  else if (parent == levels(all_parents)[6])
    return(colorRampPalette(c("hotpink", "hotpink4"))(n))
  else if (parent == levels(all_parents)[7])
    return(colorRampPalette(c("cyan", "darkcyan"))(n))
  else if (parent == levels(all_parents)[8])
    return(colorRampPalette(c("olivedrab1", "olivedrab4"))(n))
  else if (parent == levels(all_parents)[9])
    return(colorRampPalette(c("lightsalmon", "lightsalmon4"))(n))
  else if (parent == levels(all_parents)[10])
    return(colorRampPalette(c("orchid1", "purple4"))(n))
  else if (parent == levels(all_parents)[11])
    return(colorRampPalette(c("antiquewhite", "antiquewhite4"))(n))
}


optimal_layout <- function(n) {
  # Start with a square layout as close as possible
  rows <- floor(sqrt(n))
  cols <- ceiling(n / rows)

  # Adjust rows and cols if necessary
  while (rows * cols < n) {
    rows <- rows + 1
  }

  return(c(rows, cols))
}
