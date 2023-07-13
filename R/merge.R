#' Merge Multiple Genomic Datasets
#'
#' This function takes any number of genomic datasets of class 'mappoly2.data' and merges them.
#' If only one dataset is provided, the function will issue a warning and return the original dataset.
#' All datasets need to be of class 'mappoly2.data'.
#' Additional data screening and filtering is performed, including computing chi-square p-values,
#' screening non-conforming markers, and filtering redundant markers.
#'
#' @param ... Comma-separated list of datasets of class 'mappoly2.data'.
#'
#' @param filter.non.conforming if \code{TRUE} (default), data points with
#' unexpected genotypes (i.e. double reduction) are converted to 'NA'.
#' See the \code{\link[mappoly]{segreg_poly}} function for information on
#' expected classes and their respective frequencies.
#'
#' @param filter.redundant logical. If \code{TRUE} (default), removes redundant
#' markers during map construction, keeping them annotated to export to the
#' final map.
#'
#' @param verbose if \code{TRUE} (default), shows the current progress; if
#' \code{FALSE}, no output is produced
#'
#' @return The merged and filtered dataset if more than one dataset was provided;
#'         the original dataset if only one was provided.
#'
#' @examples
#' \dontrun{
#' data1 <- subset(B2721, type = "marker", n = 200)
#' data1 <- subset(data1, type = "individual", n = 80)
#' data2 <- subset(B2721, type = "marker", n = 300)
#' data2 <- subset(data2, type = "individual", n = 100)
#' data3 <- subset(B2721, type = "marker", n = 400)
#' data3 <- subset(data3, type = "individual", n = 100)
#'
#' merged_data <- merge_multiple_datasets(data1, data2, data3)
#'
#' merged_data
#'
#' plot(merged_data)
#' }
#'
#' @export
merge_multiple_datasets <- function(..., filter.non.conforming = TRUE,
                                    filter.redundant = TRUE,
                                    verbose = TRUE){

  # Put all datasets in a list
  datasets <- list(...)

  # Check if only one dataset is provided
  if(length(datasets) == 1){
    warning("Only one dataset provided. The original dataset will be returned.")
    return(datasets[[1]])
  }

  # Check that all datasets are of the correct class
  if(any(sapply(datasets, function(x) class(x) != "mappoly2.data"))){
    stop("All datasets need to be of class 'mappoly2.data'")
  }
  if (verbose) cat(" -->  Merging datasets.\n     ")

  # Use Reduce function to iteratively merge all datasets
  merged_res <- Reduce(function(x, y) merge_datasets(x, y), datasets)

  # Computing chi-square p.values
  res <- suppressWarnings(mappoly_chisq_test(merged_res))

  # Screening non-conforming markers
  if (filter.non.conforming) {
    if (verbose) cat(" -->  Filtering non-conforming markers.\n     ")
    res <- filter_non_conforming_classes(res)
  }
  # FIXME
  class(res) <- class(merged_res)
  # Screening redundant markers
  if(filter.redundant){
    if (verbose) cat(" -->  Filtering redundant markers.\n     ")
    s <- make_sequence(res, arg = "all")
    sf <- filter_redundant(s)
    res$redundant <- sf$redundant
    res <- subset_data(res, select.mrk = setdiff(res$mrk.names, sf$redundant$removed))
  }
  if(verbose) cat("----------------------------------\n")
  return(res)
}

# Function to merge two datasets
merge_data <- function(res1, res2){

  #FIXME: check for null elements

  # Ensure the class of both objects is 'mappoly2.data'
  if(class(res1) != "mappoly2.data" | class(res2) != "mappoly2.data"){
    stop("Both datasets need to be of class 'mappoly2.data'")
  }

  # Identify union of markers and individuals
  union_markers <- union(names(res1$dosage.p1), names(res2$dosage.p1))
  union_individuals <- union(res1$ind.names, res2$ind.names)

  # Check if intersecting markers have same dosage in both parents
  intersect_markers <- intersect(names(res1$dosage.p1), names(res2$dosage.p1))
  equal_dosage_p1 <- res1$dosage.p1[intersect_markers] == res2$dosage.p1[intersect_markers]
  equal_dosage_p2 <- res1$dosage.p2[intersect_markers] == res2$dosage.p2[intersect_markers]

  # Discard markers with different dosages & issue warning
  markers_to_discard <- intersect_markers[!equal_dosage_p1 | !equal_dosage_p2]
  if(length(markers_to_discard) > 0){
    warning("Discarded markers with different dosages in parents: ", markers_to_discard)
    union_markers <- setdiff(union_markers, markers_to_discard)
  }

  # Combine geno.dose and other marker-related data
  # Initialize new data objects with NA for all markers/individuals
  combined_geno_dose <- matrix(NA, nrow = length(union_markers), ncol = length(union_individuals),
                               dimnames = list(union_markers, union_individuals))
  combined_dosage_p1 <- combined_dosage_p2 <- setNames(rep(NA, length(union_markers)), union_markers)
  combined_chrom <- combined_genome_pos <- combined_ref <- combined_alt <- combined_dosage_p1

  # Fill in data from res1 and res2 where available
  combined_geno_dose[res1$mrk.names, res1$ind.names] <- res1$geno.dose
  combined_geno_dose[res2$mrk.names, res2$ind.names] <- res2$geno.dose
  combined_dosage_p1[names(res1$dosage.p1)] <- res1$dosage.p1
  combined_dosage_p1[names(res2$dosage.p1)] <- res2$dosage.p1
  combined_dosage_p2[names(res1$dosage.p2)] <- res1$dosage.p2
  combined_dosage_p2[names(res2$dosage.p2)] <- res2$dosage.p2
  combined_chrom[names(res1$chrom)] <- res1$chrom
  combined_chrom[names(res2$chrom)] <- res2$chrom
  combined_genome_pos[names(res1$genome.pos)] <- res1$genome.pos
  combined_genome_pos[names(res2$genome.pos)] <- res2$genome.pos
  combined_ref[names(res1$ref)] <- res1$ref
  combined_ref[names(res2$ref)] <- res2$ref
  combined_alt[names(res1$alt)] <- res1$alt
  combined_alt[names(res2$alt)] <- res2$alt

  # Function to extract numeric portion from chromosome names
  extract_numeric <- function(chrom) {
    as.numeric(gsub("\\D", "", chrom))
  }

  # Create ordering
  ordered_index <- order(sapply(combined_chrom, extract_numeric), combined_genome_pos)

  # Apply ordering
  combined_geno_dose <- combined_geno_dose[ordered_index,]
  combined_dosage_p1 <- combined_dosage_p1[ordered_index]
  combined_dosage_p2 <- combined_dosage_p2[ordered_index]
  combined_chrom <- combined_chrom[ordered_index]
  combined_genome_pos <- combined_genome_pos[ordered_index]
  combined_ref <- combined_ref[ordered_index]
  combined_alt <- combined_alt[ordered_index]

  # Create new merged dataset
  merged_res <- list(
    ploidy.p1 = res1$ploidy.p1,
    ploidy.p2 = res1$ploidy.p2,
    n.ind = length(union_individuals),
    n.mrk = length(union_markers),
    ind.names = union_individuals,
    mrk.names = union_markers[ordered_index], # apply ordering to marker names
    name.p1 = res1$name.p1,
    name.p2 = res1$name.p2,
    dosage.p1 = combined_dosage_p1,
    dosage.p2 = combined_dosage_p2,
    chrom = combined_chrom,
    genome.pos = combined_genome_pos,
    ref = combined_ref,
    alt = combined_alt,
    all.mrk.depth = c(res1$all.mrk.depth, res2$all.mrk.depth),
    geno.dose = combined_geno_dose,
    redundant = NULL
  )

  # Assign class to the new object
  class(merged_res) <- "mappoly2.data"

  return(merged_res)
}
