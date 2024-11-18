#' Merge Multiple Genomic Datasets
#'
#' This function takes any number of genomic datasets of class \code{mappoly2.data} and merges them.
#' If only one dataset is provided, the function will issue a warning and return the original dataset.
#' All datasets need to be of class \code{mappoly2.data}.
#' Additional data screening and filtering are performed, including computing chi-square p-values,
#' screening non-conforming markers, and filtering redundant markers.
#'
#' @param ... Comma-separated list of datasets of class \code{mappoly2.data}.
#' @param filter.non.conforming Logical. If \code{TRUE} (default), data points with
#'   unexpected genotypes (e.g., double reduction) are converted to \code{NA}.
#'   See the \code{\link[mappoly]{segreg_poly}} function for information on
#'   expected classes and their respective frequencies.
#' @param filter.redundant Logical. If \code{TRUE} (default), removes redundant
#'   markers during map construction, keeping them annotated to export to the
#'   final map.
#' @return The merged and filtered dataset if more than one dataset was provided;
#'   the original dataset if only one was provided.
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
#' merged_data <- merge_datasets(data1, data2, data3)
#'
#' merged_data
#'
#' plot(merged_data)
#' }
#' @importFrom stats setNames
#' @export
merge_datasets <- function(..., filter.non.conforming = TRUE,
                           filter.redundant = TRUE){

  # Put all datasets in a list
  datasets <- list(...)

  # Check if only one dataset is provided
  if(length(datasets) == 1){
    warning("Only one dataset provided. The original dataset will be returned.")
    return(datasets[[1]])
  }

  # Check that all datasets are of the correct class
  assert_that(all(sapply(datasets, function(x) inherits(x, "mappoly2.data"))))

  # Use Reduce function to iteratively merge all datasets
  res <- Reduce(function(x, y) merge(x, y), datasets)

  # Screening non-conforming markers
  if (filter.non.conforming) {
    res <- filter_non_conforming_classes(res)
  }

  if(filter.redundant){
    res <- suppressMessages(filter_redundant(res))
  }

  # Update QAQC values
  res$QAQC.values <- .setQAQC(id.mrk = res$mrk.names,
                              id.ind = res$ind.names,
                              miss.mrk = rowMeans(is.na(res$geno.dose)),
                              miss.ind = colMeans(is.na(res$geno.dose)),
                              chisq.pval = suppressWarnings(mappoly_chisq_test(res)))

  # Update metadata
  # Collect metrics from the datasets
  total_markers_before <- sum(sapply(datasets, function(d) d$n.mrk))
  total_individuals_before <- sum(sapply(datasets, function(d) length(d$ind.names)))
  unique_markers_before <- length(unique(unlist(sapply(datasets, function(d) d$mrk.names))))
  unique_individuals_before <- length(unique(unlist(sapply(datasets, function(d) d$ind.names))))
  total_markers_after <- res$n.mrk
  total_individuals_after <- res$n.ind

  # Prepare map_step data frame
  map_step <- data.frame(
    Metric = c("Total markers before merging",
               "Unique markers before merging",
               "Total individuals before merging",
               "Unique individuals before merging",
               "Total markers after merging",
               "Total individuals after merging",
               "Filter non-conforming markers",
               "Filter redundant markers"),
    Value = c(total_markers_before,
              unique_markers_before,
              total_individuals_before,
              unique_individuals_before,
              total_markers_after,
              total_individuals_after,
              filter.non.conforming,
              filter.redundant),
    stringsAsFactors = FALSE
  )

  # Update metadata
  res <- update_metadata(res,
                         map_step = map_step,
                         class_suffix = NULL)

  return(res)
}

#' @export
merge.mappoly2.data <- function(x, y, ...){

  # Identify union of markers and individuals
  union_markers <- union(names(x$dosage.p1), names(y$dosage.p1))
  union_individuals <- union(x$ind.names, y$ind.names)

  # Check if intersecting markers have same dosage in both parents
  intersect_markers <- intersect(names(x$dosage.p1), names(y$dosage.p1))
  equal_dosage_p1 <- x$dosage.p1[intersect_markers] == y$dosage.p1[intersect_markers]
  equal_dosage_p2 <- x$dosage.p2[intersect_markers] == y$dosage.p2[intersect_markers]

  # Discard markers with different dosages & issue warning
  markers_to_discard <- intersect_markers[!equal_dosage_p1 | !equal_dosage_p2]
  if(length(markers_to_discard) > 0){
    warning("Discarded markers with different dosages in parents: ", paste(markers_to_discard, collapse = ", "))
    union_markers <- setdiff(union_markers, markers_to_discard)
  }

  # Combine geno.dose and other marker-related data
  # Initialize new data objects with NA for all markers/individuals
  combined_geno_dose <- matrix(NA, nrow = length(union_markers), ncol = length(union_individuals),
                               dimnames = list(union_markers, union_individuals))
  combined_dosage_p1 <- combined_dosage_p2 <- setNames(rep(NA, length(union_markers)), union_markers)
  combined_chrom <- combined_genome_pos <- combined_ref <- combined_alt <- combined_dosage_p1

  # Fill in data from x and y where available
  combined_geno_dose[x$mrk.names, x$ind.names] <- x$geno.dose
  combined_geno_dose[y$mrk.names, y$ind.names] <- y$geno.dose
  combined_dosage_p1[names(x$dosage.p1)] <- x$dosage.p1
  combined_dosage_p1[names(y$dosage.p1)] <- y$dosage.p1
  combined_dosage_p2[names(x$dosage.p2)] <- x$dosage.p2
  combined_dosage_p2[names(y$dosage.p2)] <- y$dosage.p2
  combined_chrom[names(x$chrom)] <- x$chrom
  combined_chrom[names(y$chrom)] <- y$chrom
  combined_genome_pos[names(x$genome.pos)] <- x$genome.pos
  combined_genome_pos[names(y$genome.pos)] <- y$genome.pos
  combined_ref[names(x$ref)] <- x$ref
  combined_ref[names(y$ref)] <- y$ref
  combined_alt[names(x$alt)] <- x$alt
  combined_alt[names(y$alt)] <- y$alt

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
    ploidy.p1 = x$ploidy.p1,
    ploidy.p2 = x$ploidy.p2,
    n.ind = length(union_individuals),
    n.mrk = length(union_markers),
    ind.names = union_individuals,
    mrk.names = union_markers[ordered_index], # apply ordering to marker names
    name.p1 = x$name.p1,
    name.p2 = x$name.p2,
    dosage.p1 = combined_dosage_p1,
    dosage.p2 = combined_dosage_p2,
    chrom = combined_chrom,
    genome.pos = combined_genome_pos,
    ref = combined_ref,
    alt = combined_alt,
    all.mrk.depth = c(x$all.mrk.depth, y$all.mrk.depth),
    geno.dose = combined_geno_dose,
    redundant = NULL,
    QAQC.values = NULL  # Will be updated later
  )

  # Assign class to the new object
  class(merged_res) <- "mappoly2.data"

  # Since the merge is an internal step, we'll let the main function update the metadata
  # So we won't call update_metadata here

  return(merged_res)
}
