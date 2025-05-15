#' Read Genetic Marker Data from a CSV File
#'
#' Reads genetic marker data from a CSV file and returns an object of class \code{mappoly2.data}.
#'
#' The CSV file should have markers in rows, with the first row serving as the header.
#' The first seven columns must contain the following information, in order:
#' \enumerate{
#'   \item Marker names
#'   \item Dosage in parent 1
#'   \item Dosage in parent 2
#'   \item Chromosome information (e.g., chromosome, scaffold, contig)
#'   \item Position of the marker within the sequence
#'   \item Alternate allele (if available)
#'   \item Reference allele (if available)
#' }
#' If allele information is not available, the values should be set to \code{NA}.
#' The remaining columns should contain the dosage for each individual in the full-sib population.
#' Refer to the \code{Examples} section for a tetraploid example.
#'
#' @param file.in A character string specifying the name or full path to the input CSV file.
#' @param ploidy.p1 The ploidy level of parent 1.
#' @param ploidy.p2 The ploidy level of parent 2 (defaults to \code{ploidy.p1}).
#' @param name.p1 The name of parent 1.
#' @param name.p2 The name of parent 2.
#' @param filter.non.conforming Logical value indicating whether to convert non-conforming data points
#'   (e.g., double reduction) to \code{NA} (default is \code{TRUE}).
#' @param filter.redundant Logical value indicating whether to remove redundant markers during map construction,
#'   retaining annotations for export in the final map (default is \code{TRUE}).
#' @param verbose Logical value indicating whether to display progress updates (default is \code{TRUE}).
#'
#' @return An object of class \code{mappoly2.data}, which is a list containing the following components:
#' \describe{
#'   \item{ploidy.p1}{Ploidy level of parent 1.}
#'   \item{ploidy.p2}{Ploidy level of parent 2.}
#'   \item{n.ind}{Number of individuals in the population.}
#'   \item{n.mrk}{Total number of markers.}
#'   \item{ind.names}{Names or identifiers of the individuals.}
#'   \item{mrk.names}{Names or identifiers of the genetic markers.}
#'   \item{name.p1}{Name or identifier of parent 1.}
#'   \item{name.p2}{Name or identifier of parent 2.}
#'   \item{dosage.p1}{Dosage for parent 1.}
#'   \item{dosage.p2}{Dosage for parent 2.}
#'   \item{chrom}{Chromosome numbers for all markers.}
#'   \item{genome.pos}{Physical positions of the markers within the genome.}
#'   \item{ref}{Reference alleles for the markers.}
#'   \item{alt}{Alternate alleles for the markers.}
#'   \item{all.mrk.depth}{Depth of coverage for all markers (NULL when using CSV input files).}
#'   \item{geno.dose}{Matrix of dosages for each marker (rows) and each individual (columns).}
#'   \item{redundant}{List of non-redundant markers and their equivalent redundant markers if
#'     \code{filter.redundant} is \code{TRUE}.}
#'   \item{QAQC.values}{List containing quality assurance and quality control values:
#'     \describe{
#'       \item{markers}{Data frame with statistics for each marker, including \code{miss}
#'         (missing data rate) and \code{chisq.pval} (chi-squared test p-value).}
#'       \item{individuals}{Data frame with statistics for each individual, including \code{miss}
#'         (missing data rate) and \code{full.sib} (indicator of non-belonging to the analyzed
#'         bi-parental cross).}
#'     }
#'   }
#' }
#'
#' @examples
#' \donttest{
#' # Read a tetraploid dataset from a CSV file
#' tempfl <- list.files(system.file('extdata', package = 'mappoly2'), full.names = TRUE)
#' alfalfa.bc <- read_geno_csv(
#'   file.in = tempfl,
#'   ploidy.p1 = 4,
#'   name.p1 = "I195",
#'   name.p2 = "F1.85.209"
#' )
#' print(alfalfa.bc, detailed = TRUE)
#' plot(alfalfa.bc)
#' }
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @importFrom grDevices rgb
#' @importFrom stats na.omit
#' @importFrom utils read.csv
#' @importFrom assertthat is.readable
#' @export
read_geno_csv <- function(
    file.in,
    ploidy.p1,
    ploidy.p2 = ploidy.p1,
    name.p1 = NULL,
    name.p2 = NULL,
    filter.non.conforming = TRUE,
    filter.redundant = TRUE,
    verbose = TRUE
) {
  # Ensure the input file is readable
  assert_that(is.readable(file.in))

  # Read the CSV file
  data <- read.csv(file = file.in, header = TRUE, stringsAsFactors = FALSE)

  # Store the number of original markers and individuals
  num_markers <- nrow(data)
  num_individuals <- ncol(data) - 7  # Adjusted for the first seven columns

  # Convert the data to mappoly format
  mappoly_data <- mappoly2:::table_to_mappoly(
    data,
    ploidy.p1 = ploidy.p1,
    ploidy.p2 = ploidy.p2,
    name.p1 = name.p1,
    name.p2 = name.p2,
    filter.non.conforming = filter.non.conforming,
    filter.redundant = filter.redundant,
    verbose = verbose
  )

  # Define metrics for metadata
  metrics <- c(
    "Input file",
    "MD5 hash",
    "File size (MB)",
    "Date",
    "Original markers from file",
    "Original individuals from file"
  )

  # Define the corresponding values
  values <- c(
    file.in,
    tools::md5sum(file.in),
    round(file.size(file.in) / 1024000, 2),
    date(),
    num_markers,
    num_individuals
  )

  # Create the map_step data frame
  map_step <- data.frame(
    Metric = metrics,
    Value = values,
    stringsAsFactors = FALSE
  )

  # Update the metadata of the mappoly object
  mappoly_data <- update_metadata(
    mappoly_data,
    map_step = map_step,
    class_suffix = "mappoly.init.filter"
  )

  return(mappoly_data)
}

#' Read Genetic Marker Data from a VCF File
#'
#' Reads genetic marker data from a Variant Call Format (VCF) file and returns an object of class
#' \code{mappoly2.data}.
#'
#' Supports VCF files of version 4.0 or higher.
#'
#' @param file.in A character string specifying the name or full path of the input VCF file.
#' @param ploidy.p1 The ploidy level of parent 1.
#' @param ploidy.p2 The ploidy level of parent 2 (defaults to \code{ploidy.p1}).
#' @param name.p1 The name of parent 1.
#' @param name.p2 The name of parent 2.
#' @param name.offspring A character vector containing the names of the offspring (defaults to all
#'   individuals except the parents).
#' @param filter.non.conforming Logical value indicating whether to convert non-conforming data points
#'   to \code{NA} (default is \code{TRUE}).
#' @param filter.redundant Logical value indicating whether to remove redundant markers during map
#'   construction, retaining annotations for export in the final map (default is \code{TRUE}).
#' @param verbose Logical value indicating whether to display progress information (default is \code{TRUE}).
#' @param min.gt.depth Minimum genotype depth to retain information; values below this threshold are
#'   replaced with \code{NA} (default is 0).
#' @param min.av.depth Minimum average depth to retain markers (default is 0).
#' @param max.missing Maximum proportion of missing data allowed to retain markers (range from 0 to 1;
#'   default is 1).
#'
#' @return An object of class \code{mappoly2.data}, which is a list containing the following components:
#' \describe{
#'   \item{ploidy.p1}{Ploidy level of parent 1.}
#'   \item{ploidy.p2}{Ploidy level of parent 2.}
#'   \item{n.ind}{Number of individuals.}
#'   \item{n.mrk}{Total number of markers.}
#'   \item{ind.names}{Names of the individuals.}
#'   \item{mrk.names}{Names of the genetic markers.}
#'   \item{name.p1}{Name of parent 1.}
#'   \item{name.p2}{Name of parent 2.}
#'   \item{dosage.p1}{Dosage information for parent 1.}
#'   \item{dosage.p2}{Dosage information for parent 2.}
#'   \item{chrom}{Chromosome numbers for all markers.}
#'   \item{genome.pos}{Physical positions of the markers within the genome.}
#'   \item{ref}{Reference alleles for the markers.}
#'   \item{alt}{Alternate alleles for the markers.}
#'   \item{all.mrk.depth}{Depth of coverage for all markers.}
#'   \item{geno.dose}{Matrix of dosages for each marker (rows) and each individual (columns).}
#'   \item{redundant}{List of non-redundant markers and their equivalent redundant markers if
#'     \code{filter.redundant} is \code{TRUE}.}
#'   \item{QAQC.values}{Quality assurance and quality control values:
#'     \describe{
#'       \item{markers}{Data frame with statistics for each marker, including \code{miss}
#'         (missing data rate) and \code{chisq.pval} (chi-squared test p-value).}
#'       \item{individuals}{Data frame with statistics for each individual, including \code{miss}
#'         (missing data rate) and \code{full.sib} (indicator of non-belonging to the analyzed
#'         bi-parental cross).}
#'     }
#'   }
#' }
#'
#' @examples
#' \donttest{
#' # Hexaploid sweetpotato: Subset of chromosome 3
#' fl <- "https://github.com/mmollina/MAPpoly_vignettes/raw/master/data/sweet_sample_ch3.vcf.gz"
#' tempfl <- tempfile(pattern = 'chr3_', fileext = '.vcf.gz')
#' download.file(fl, destfile = tempfl)
#' dat.dose.vcf <- read_vcf(
#'   file.in = tempfl,
#'   ploidy.p1 = 6,
#'   name.p1 = "PARENT1",
#'   name.p2 = "PARENT2",
#'   min.av.depth = 20
#' )
#' print(dat.dose.vcf)
#' plot(dat.dose.vcf)
#' }
#'
#' @references Gabriel Gesteira, \email{gdesiqu@ncsu.edu}; Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
read_vcf <- function(
    file.in,
    ploidy.p1,
    ploidy.p2 = ploidy.p1,
    name.p1,
    name.p2,
    name.offspring = NULL,
    filter.non.conforming = TRUE,
    filter.redundant = TRUE,
    verbose = TRUE,
    min.gt.depth = 0,
    min.av.depth = 0,
    max.missing = 1
) {
  # Ensure the input file is readable
  assert_that(is.readable(file.in))

  # Validate ploidy levels
  mappoly2:::has.adequate.ploidy(ploidy.p1)
  mappoly2:::has.adequate.ploidy(ploidy.p2)

  # Expected ploidy level for offspring
  expected_ploidy_offspring <- (ploidy.p1 + ploidy.p2) / 2

  # Get the full path and size of the input file
  input_file <- normalizePath(file.in)
  input_size_mb <- file.size(file.in) / 1024^2  # File size in MB

  # Warn if the file size is too large
  if (verbose && input_size_mb > 3000) {
    warning("Your VCF file is larger than 3 GB. Check for available RAM memory.")
  }

  if (verbose) cat("Reading data...\n")

  # Read the VCF file
  vcf_data <- vcfR::read.vcfR(file.in, verbose = verbose)

  # Get individual names from the genotype matrix
  ind_names <- colnames(vcf_data@gt)[-1]

  # Check if parent names are in the dataset
  if (!name.p1 %in% ind_names) stop(name.p1, " is not in the dataset")
  if (!name.p2 %in% ind_names) stop(name.p2, " is not in the dataset")

  # Determine offspring names
  if (is.null(name.offspring)) name.offspring <- setdiff(ind_names, c(name.p1, name.p2))
  name.offspring <- intersect(name.offspring, ind_names)
  if (length(name.offspring) == 0) stop("No offspring individuals found in the dataset")

  # Get number of markers and individuals
  num_markers <- nrow(vcf_data@gt)
  num_individuals <- length(name.offspring)

  # Get marker names and chromosome information
  marker_names <- vcf_data@fix[, "ID"]
  chrom <- vcf_data@fix[, "CHROM"]
  genome_pos <- as.numeric(vcf_data@fix[, "POS"])

  # Handle unnamed markers
  unnamed_markers <- is.na(marker_names) | marker_names == ""
  if (any(unnamed_markers)) {
    if (verbose) {
      cat("No marker names found. Generating names using chromosome and position information.\n")
    }
    marker_names[unnamed_markers] <- paste0(chrom[unnamed_markers], "_", genome_pos[unnamed_markers])
  }

  # Assign names to chromosome and position vectors
  names(chrom) <- marker_names
  names(genome_pos) <- marker_names

  # Get reference and alternate alleles
  ref_alleles <- vcf_data@fix[, "REF"]
  alt_alleles <- vcf_data@fix[, "ALT"]
  names(ref_alleles) <- names(alt_alleles) <- marker_names

  if (verbose) cat("Processing genotypes...\n")

  # Identify positions of GT (genotype) and DP (depth) fields
  format_fields <- unlist(strsplit(vcf_data@gt[1, 1], ":"))
  gt_pos <- which(format_fields == "GT")
  dp_pos <- which(format_fields == "DP")

  # Get ploidy levels and depths for all individuals
  offspring_ploidy <- mappoly2:::.vcf_get_ploidy(vcf_data@gt[, name.offspring], gt_pos)
  geno_depth <- mappoly2:::.vcf_get_depth(
    vcf_data@gt[, c(name.p1, name.p2, name.offspring)],
    dp_pos
  )

  # Assign dimnames to geno_depth
  dimnames(geno_depth) <- list(marker_names, c(name.p1, name.p2, name.offspring))

  # Transform genotype data to dosage information
  geno_dose <- mappoly2:::.vcf_transform_dosage(
    vcf_data@gt[, c(name.p1, name.p2, name.offspring)],
    gt_pos
  )
  dimnames(geno_dose) <- list(marker_names, c(name.p1, name.p2, name.offspring))
  geno_dose[geno_dose == -1] <- NA

  if (verbose) cat("Done processing genotypes.\n")

  # Check if observed ploidy matches expected ploidy
  unique_offspring_ploidy <- unique(c(offspring_ploidy))
  if (!expected_ploidy_offspring %in% unique_offspring_ploidy) {
    stop("Observed ploidy in offspring differs from expected ploidy.")
  }

  # Filter markers based on ploidy, depth, and missing data thresholds
  valid_ploidy_markers <- rowSums(offspring_ploidy == expected_ploidy_offspring) == num_individuals
  avg_marker_depth <- rowMeans(geno_depth)
  sufficient_depth_markers <- avg_marker_depth >= min.av.depth
  acceptable_missing_markers <- rowSums(is.na(geno_dose)) / ncol(geno_dose) <= max.missing

  selected_markers <- marker_names[
    valid_ploidy_markers & sufficient_depth_markers & acceptable_missing_markers
  ]

  # Update genotype dosage and depth matrices
  geno_dose <- geno_dose[selected_markers, , drop = FALSE]
  geno_depth <- geno_depth[selected_markers, , drop = FALSE]
  geno_dose[geno_depth < min.gt.depth] <- NA
  avg_marker_depth <- avg_marker_depth[selected_markers]

  # Ensure parental genotypes are available
  parental_genotypes_complete <- rowSums(
    is.na(geno_dose[, c(name.p1, name.p2)])
  ) == 0
  markers_with_parents <- rownames(geno_dose)[parental_genotypes_complete]

  # Create a data frame for the selected markers
  marker_data <- data.frame(
    snp_id = markers_with_parents,
    P1 = as.integer(geno_dose[markers_with_parents, name.p1]),
    P2 = as.integer(geno_dose[markers_with_parents, name.p2]),
    chrom = chrom[markers_with_parents],
    genome_pos = genome_pos[markers_with_parents],
    ref = ref_alleles[markers_with_parents],
    alt = alt_alleles[markers_with_parents],
    geno_dose[markers_with_parents, name.offspring, drop = FALSE],
    stringsAsFactors = FALSE
  )

  # Convert the data frame to a mappoly object
  mappoly_data <- mappoly2:::table_to_mappoly(
    marker_data,
    ploidy.p1 = ploidy.p1,
    ploidy.p2 = ploidy.p2,
    name.p1 = name.p1,
    name.p2 = name.p2,
    filter.non.conforming = filter.non.conforming,
    filter.redundant = filter.redundant,
    verbose = verbose
  )

  # Add read depth to the QAQC values
  mappoly_data$QAQC.values$markers$read.depth <- avg_marker_depth[
    rownames(mappoly_data$QAQC.values$markers)
  ]

  # Update metadata
  metrics <- c(
    "Input file",
    "MD5 hash",
    "File size (MB)",
    "Date",
    "Original markers from file",
    "Original individuals from file",
    "Markers with expected offspring ploidy",
    "Minimum average depth",
    "Markers passing depth filter",
    "Maximum missing threshold",
    "Markers passing missing filter",
    "Selected markers after screening"
  )

  values <- c(
    file.in,
    tools::md5sum(file.in),
    round(file.size(file.in) / 1024^2, 2),
    date(),
    num_markers,
    num_individuals,
    sum(valid_ploidy_markers),
    min.av.depth,
    sum(sufficient_depth_markers),
    max.missing,
    sum(acceptable_missing_markers),
    length(markers_with_parents)
  )

  map_step <- data.frame(
    Metric = metrics,
    Value = values,
    stringsAsFactors = FALSE
  )

  # Update the metadata of the mappoly object
  mappoly_data <- update_metadata(
    mappoly_data,
    map_step = map_step,
    class_suffix = "mappoly.init.filter"
  )

  return(mappoly_data)
}

#' Create a List of `mappoly.data` Objects for F1 Crosses
#'
#' This function takes genotype and pedigree files to generate
#' a list of `mappoly.data` objects for multiple F1 populations. It performs
#' optional filtering for non-conforming dosage classes and redundant markers.
#'
#' @param data.file A string. Path to the CSV file with genotype dosage data, including marker ID,
#'   parental genotypes, marker metadata (e.g., chromosome, position, ref/alt),
#'   and progeny dosage calls.
#' @param pedigree.file A string. Path to the CSV file with three columns: `cross`, `parent_1`, and `parent_2`,
#'   indicating the cross name and its respective parents.
#' @param ploidy.vec A named numeric vector indicating the ploidy level for each parent.
#'   Names must correspond to parental names in the `pedigree.file`.
#' @param filter.non.conforming Logical. If `TRUE`, remove genotype calls that violate
#'   the expected segregation range under no double reduction.
#' @param filter.redundant Logical. If `TRUE`, remove redundant markers.
#' @param verbose Logical. If `TRUE`, print progress and summary messages.
#'
#' @return An object of class `mappoly.data.list`, a named list of `mappoly.data` objects
#'   for each cross found in the pedigree.
#'
#' @export
#'
#' @examples
#' ploidy.vec <- c(4, 2, 4, 2, 4, 4)
#' names(ploidy.vec) <- c("P1", "P2", "P3", "P4", "P5", "P6")
#' result <- read_multipop_geno_csv(
#'   data.file = "path/to/data_multi_pop.csv",
#'   pedigree.file = "path/to/pedigree.csv",
#'   ploidy.vec = ploidy.vec
#' )
read_multipop_geno_csv <- function(data.file,
                                     pedigree.file,
                                     ploidy.vec,
                                     filter.non.conforming = TRUE,
                                     filter.redundant = TRUE,
                                     verbose = TRUE) {
  dat <- readr::read_csv(data.file, show_col_types = FALSE)
  pedigree <- readr::read_csv(pedigree.file, show_col_types = FALSE)

  dat_split <- mappoly2:::create_all_dat(dat, pedigree)

  dat_split <- lapply(dat_split, function(df) {
    df %>%
      dplyr::select(
        1,        # snp_id
        6, 7,     # parent columns (P1 and P2)
        2:5,      # chrom, genome_pos, ref, alt
        8:ncol(.) # progeny genotypes
      )
  })

  crosses <- names(dat_split)
  out_list <- vector("list", length(crosses))
  names(out_list) <- crosses

  for (i in crosses) {
    cur_f1 <- dplyr::filter(pedigree, cross == i)
    p1 <- unique(cur_f1$parent_1)
    p2 <- unique(cur_f1$parent_2)
    stopifnot(length(p1) == 1, length(p2) == 1)

    out_list[[i]] <- mappoly2:::table_to_mappoly(
      x = as.data.frame(dat_split[[i]]),
      ploidy.p1 = ploidy.vec[p1],
      ploidy.p2 = ploidy.vec[p2],
      name.p1 = p1,
      name.p2 = p2,
      filter.non.conforming = filter.non.conforming,
      filter.redundant = filter.redundant,
      verbose = verbose
    )

    if (verbose) {
      message("Processed cross: ", i)
    }
  }

  structure(out_list, class = "mappoly.data.list")
}

#' Create a List of `mappoly.data` Objects for F1 Crosses
#'
#' This function takes genotype and pedigree files to generate
#' a list of `mappoly.data` objects for multiple F1 populations. It performs
#' optional filtering for non-conforming dosage classes and redundant markers.
#'
#' @param data.file A string. Path to the CSV file with genotype dosage data, including marker ID,
#'   parental genotypes, marker metadata (e.g., chromosome, position, ref/alt),
#'   and progeny dosage calls.
#' @param pedigree.file A string. Path to the CSV file with three columns: `cross`, `parent_1`, and `parent_2`,
#'   indicating the cross name and its respective parents.
#' @param ploidy.vec A named numeric vector indicating the ploidy level for each parent.
#'   Names must correspond to parental names in the `pedigree.file`.
#' @param filter.non.conforming Logical. If `TRUE`, remove genotype calls that violate
#'   the expected segregation range under no double reduction.
#' @param filter.redundant Logical. If `TRUE`, remove redundant markers.
#' @param verbose Logical. If `TRUE`, print progress and summary messages.
#'
#' @return An object of class `mappoly.data.list`, a named list of `mappoly.data` objects
#'   for each cross found in the pedigree.
#'
#' @export
#' @importFrom readr read_csv
#' @examples
#' ploidy.vec <- c(4, 2, 4, 2, 4, 4)
#' names(ploidy.vec) <- c("P1", "P2", "P3", "P4", "P5", "P6")
#' result <- read_multipop_geno_csv(
#'   data.file = "path/to/data_multi_pop.csv",
#'   pedigree.file = "path/to/pedigree.csv",
#'   ploidy.vec = ploidy.vec
#' )
read_multipop_geno_csv <- function(data.file,
                                     pedigree.file,
                                     ploidy.vec,
                                     filter.non.conforming = TRUE,
                                     filter.redundant = TRUE,
                                     verbose = TRUE) {
  dat <- read_csv(data.file, show_col_types = FALSE)
  pedigree <- read_csv(pedigree.file, show_col_types = FALSE)

  dat_split <- create_all_dat(dat, pedigree)

  dat_split <- lapply(dat_split, function(df) {
    df %>%
      dplyr::select(
        1,        # snp_id
        6, 7,     # parent columns (P1 and P2)
        2:5,      # chrom, genome_pos, ref, alt
        8:ncol(.) # progeny genotypes
      )
  })

  crosses <- names(dat_split)
  out_list <- vector("list", length(crosses))
  names(out_list) <- crosses

  for (i in crosses) {
    cur_f1 <- dplyr::filter(pedigree, cross == i)
    p1 <- unique(cur_f1$parent_1)
    p2 <- unique(cur_f1$parent_2)
    stopifnot(length(p1) == 1, length(p2) == 1)

    out_list[[i]] <- mappoly2:::table_to_mappoly(
      x = as.data.frame(dat_split[[i]]),
      ploidy.p1 = ploidy.vec[p1],
      ploidy.p2 = ploidy.vec[p2],
      name.p1 = p1,
      name.p2 = p2,
      filter.non.conforming = filter.non.conforming,
      filter.redundant = filter.redundant,
      verbose = verbose
    )

    if (verbose) {
      message("Processed cross: ", i)
    }
  }

  structure(out_list, class = "mappoly2.data.list")
}


#' @export
print.mappoly2.data.list  <- function(x,...){
  y <- table(sapply(x, function(x) x$name.p1),
             sapply(x, function(x) x$name.p2))
  for(i in seq_along(x))
    y[x[[i]]$name.p1, x[[i]]$name.p2] <- x[[i]]$n.ind

  sets <- lapply(x, function(x) x$mrk.names)
  n.unique.mrk <- length(Reduce(union, sets))
  sets <- lapply(x, function(x) x$mrk.names)
  total_elements <- sum(sapply(sets, length))
  redundancy_score <- (total_elements - n.unique.mrk) / total_elements





  a1<-sapply(x$consensus.map, function(x) sapply(x$ph$PH, function(x) ncol(x)))
  a2<-sapply(x$consensus.map, function(x) is.null(x$haploprob))
  msg("Multiparental data:", col = "blue")
  cat("Ploidy of founders:             ", a1[,1], "\n")
  cat("Total No. individuals:          ", sum(y), "\n")
  cat("Total No. markers               ", sum(n.mrk), "\n")
  cat("Haplotype probability computed: ", ifelse(all(a2), "No", "Yes"), "\n\n")
  cat("Number of individuals per cross:\n")
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



#' Convert a Data Frame to mappoly2.data Object
#'
#' Converts a data frame containing genetic marker data into an object of class \code{mappoly2.data}.
#'
#' @param x A data frame containing genetic marker data.
#' @param ploidy.p1 Ploidy level of parent 1.
#' @param ploidy.p2 Ploidy level of parent 2.
#' @param name.p1 Name of parent 1 (defaults to the name in the data frame).
#' @param name.p2 Name of parent 2 (defaults to the name in the data frame).
#' @param filter.non.conforming Logical value indicating whether to filter out non-conforming markers
#'   (default is \code{TRUE}).
#' @param filter.redundant Logical value indicating whether to filter out redundant markers (default is
#'   \code{TRUE}).
#' @param verbose Logical value indicating whether to display progress information (default is \code{TRUE}).
#' @keywords internal
table_to_mappoly <- function(
    x,
    ploidy.p1,
    ploidy.p2,
    name.p1 = NULL,
    name.p2 = NULL,
    filter.non.conforming = TRUE,
    filter.redundant = TRUE,
    verbose = TRUE
) {
  # Remove markers with missing parental dosages
  x <- x[!is.na(x[, 2]) & !is.na(x[, 3]), ]

  # Get the number of individuals and markers
  n.ind <- ncol(x) - 7
  n.mrk <- nrow(x)
  n.mrk.nona.on.parents <- n.mrk

  # Get marker names
  mrk.names <- as.character(x[, 1])

  # Get individual names
  ind.names <- colnames(x)[-(1:7)]

  # Get parent names
  if (is.null(name.p1)) name.p1 <- colnames(x)[2]
  if (is.null(name.p2)) name.p2 <- colnames(x)[3]

  # Get dosages for parents
  dosage.p1 <- as.integer(x[, 2])
  dosage.p2 <- as.integer(x[, 3])

  # Compute dosage differences for polymorphic markers
  d.p1 <- abs(abs(dosage.p1 - (ploidy.p1 / 2)) - (ploidy.p1 / 2))
  d.p2 <- abs(abs(dosage.p2 - (ploidy.p2 / 2)) - (ploidy.p2 / 2))
  id <- (d.p1 + d.p2) != 0
  n.segreg.mrk <- sum(id)

  # Markers with parental dosages less than or equal to ploidy levels
  within_ploidy <- (dosage.p1 <= ploidy.p1) & (dosage.p2 <= ploidy.p2)
  id <- id & within_ploidy
  n.mrk.less.equal.ploidy <- sum(id)

  # Extract chromosome and position information
  chrom <- as.character(x[, 4])
  chrom[is.na(chrom) | chrom == ""] <- "NoChr"

  genome.pos <- as.numeric(x[, 5])

  # Get reference and alternate alleles
  ref <- as.character(x[, 6])
  alt <- as.character(x[, 7])

  # Assign names to vectors
  names(ref) <- names(alt) <- names(genome.pos) <- names(chrom) <-
    names(dosage.p2) <- names(dosage.p1) <- mrk.names

  if (verbose) {
    cat("Reading data:\n")
    cat("    Ploidy level of", name.p1, ":", ploidy.p1, "\n")
    cat("    Ploidy level of", name.p2, ":", ploidy.p2, "\n")
    cat("    Number of individuals:", n.ind, "\n")
    cat("    Number of markers:", n.mrk, "\n")
    cat("    Number of informative markers:", sum(id), "(",
        round(100 * sum(id) / n.mrk, 1), "%)\n")
    if (any(!is.na(chrom))) cat("    Chromosome information is available.\n")
    if (any(!is.na(genome.pos))) cat("    Position information is available.\n")
  }

  # Get genotypic data
  geno.dose <- as.matrix(x[, -(1:7)])
  dimnames(geno.dose) <- list(mrk.names, ind.names)

  # Replace dosages exceeding offspring ploidy level with NA
  max_offspring_ploidy <- (ploidy.p1 + ploidy.p2) / 2
  geno.dose[geno.dose > max_offspring_ploidy] <- NA

  # Filter to valid markers
  geno.dose <- geno.dose[id, , drop = FALSE]

  # Create the mappoly2.data object
  mappoly_data <- list(
    ploidy.p1 = ploidy.p1,
    ploidy.p2 = ploidy.p2,
    n.ind = n.ind,
    n.mrk = sum(id),
    ind.names = ind.names,
    mrk.names = mrk.names[id],
    name.p1 = name.p1,
    name.p2 = name.p2,
    dosage.p1 = dosage.p1[id],
    dosage.p2 = dosage.p2[id],
    chrom = chrom[id],
    genome.pos = genome.pos[id],
    ref = ref[id],
    alt = alt[id],
    all.mrk.depth = NULL,
    geno.dose = geno.dose,
    redundant = NULL,
    QAQC.values = NULL
  )
  class(mappoly_data) <- c("mappoly2.data")

  # Filter non-conforming markers if requested
  if (filter.non.conforming) {
    if (verbose) cat("Filtering non-conforming markers...\n")
    mappoly_data <- mappoly2:::filter_non_conforming_classes(mappoly_data)
  }

  # Filter redundant markers if requested
  if (filter.redundant) {
    if (verbose) cat("Filtering redundant markers...\n")
    mappoly_data <- suppressMessages(filter_redundant(mappoly_data, plot = FALSE))
  }

  # Compute QAQC values
  mappoly_data$QAQC.values <- .setQAQC(
    id.mrk = mappoly_data$mrk.names,
    id.ind = mappoly_data$ind.names,
    miss.mrk = rowMeans(is.na(mappoly_data$geno.dose)),
    miss.ind = colMeans(is.na(mappoly_data$geno.dose)),
    chisq.pval = suppressWarnings(mappoly_chisq_test(mappoly_data))
  )

  if (verbose) cat("Data processing complete.\n")

  # Update metadata
  map_step <- data.frame(
    Metric = c(
      "Markers with no missing data on parents",
      "Segregating markers",
      "Markers within ploidy levels",
      paste("Ploidy level of", mappoly_data$name.p1),
      paste("Ploidy level of", mappoly_data$name.p2),
      "Number of individuals",
      "Number of markers",
      "Percentage of redundant markers",
      "Percentage of missing data",
      "Chromosome information",
      "Position information"
    ),
    Value = c(
      n.mrk.nona.on.parents,
      n.segreg.mrk,
      n.mrk.less.equal.ploidy,
      mappoly_data$ploidy.p1,
      mappoly_data$ploidy.p2,
      mappoly_data$n.ind,
      length(mappoly_data$mrk.names),
      if (!is.null(mappoly_data$redundant)) {
        paste(
          round(
            100 * nrow(mappoly_data$redundant) / (
              length(mappoly_data$mrk.names) + nrow(mappoly_data$redundant)
            ),
            1
          ),
          "%"
        )
      } else {
        "Unavailable"
      },
      paste(
        round(100 * sum(is.na(mappoly_data$geno.dose)) / length(mappoly_data$geno.dose), 1),
        "%"
      ),
      if (any(!is.na(mappoly_data$chrom))) "Available" else "Unavailable",
      if (any(!is.na(mappoly_data$genome.pos))) "Available" else "Unavailable"
    ),
    stringsAsFactors = FALSE
  )

  # Update the metadata of the mappoly object
  mappoly_data <- update_metadata(
    mappoly_data,
    map_step = map_step,
    class_suffix = NULL
  )

  return(mappoly_data)
}
