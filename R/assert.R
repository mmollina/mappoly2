#' @title Validate SNP Genotype CSV File Structure
#'
#' @description
#' Validates that a SNP genotype file (loaded as a data.frame or tibble) contains the expected structure,
#' including essential columns (e.g., SNP ID, chromosome, position, ref/alt alleles),
#' ploidy-conforming dosage data, and numeric genotype values.
#' Column names are matched flexibly using common aliases.
#'
#' @param df A data.frame or tibble containing SNP marker data.
#' @param ploidy.p1 Integer. Ploidy level of parent 1.
#' @param ploidy.p2 Integer. Ploidy level of parent 2.
#' @param max_missing_rate Numeric. Maximum allowed proportion of missing genotype data per SNP before issuing a warning. Default is `Inf`.
#' @param reorder Logical. If TRUE (default), reorders the data frame by chromosome and position if needed.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{valid}{Logical. TRUE if the dataset passed structural checks (even if warnings exist).}
#'   \item{n_snps}{Number of SNP markers (rows).}
#'   \item{n_samples}{Number of genotype/sample columns.}
#'   \item{ploidy.p1}{Ploidy of parent 1.}
#'   \item{ploidy.p2}{Ploidy of parent 2.}
#'   \item{max_na_rate_observed}{Maximum proportion of missing genotypes per SNP.}
#'   \item{problems}{Named list of warning messages issued during validation.}
#' }
#'
#' @importFrom dplyr group_by summarise filter n_distinct %>%
#' @importFrom rlang .data
#' @export
assert_valid_snp_csv <- function(df,
                                 ploidy.p1,
                                 ploidy.p2,
                                 name.p1 = NULL,
                                 name.p2 = NULL,
                                 max_missing_rate = Inf,
                                 reorder = TRUE) {
  problems <- list()

  # ---- 0. Flexible column matching setup ----
  expected_roles <- list(
    snp_id = c("snp", "snp_id", "marker", "marker_id", "snpid", "id", "locus"),
    chrom = c("chrom", "chr", "chromosome", "CHR", "CHROM", "seqid", "scaffold"),
    genome_pos = c("pos", "position", "bp", "genome_pos", "base_pair", "bp_position", "start"),
    ref = c("ref", "reference", "ref_allele", "reference_allele", "allele1"),
    alt = c("alt", "alternate", "alt_allele", "alternative", "allele2"),
    P1 = unique(tolower(name.p1)),
    P2 = unique(tolower(name.p2))
  )

  matched_cols <- list()
  cte<-1
  for (role in names(expected_roles)) {
    pattern <- expected_roles[[role]]
    match <- which(tolower(colnames(df)) %in% pattern)

    if (length(match) == 0) {
      if(role == "P1")
      {
        warning_msg <- paste("Name of parent 1 not provided. Using column", cte)
      } else if(role == "P2"){
        warning_msg <- paste("Name of parent 2 not provided. Using column", cte)
      } else {
        warning_msg <- paste("Could not find column for", role, "- expected one of:", paste(pattern, collapse = ", "),
                             " Using column", cte)
      }
      warning(warning_msg)
      problems[[role]] <- warning_msg
      matched_cols[[role]] <- NA
    } else if (length(match) > 1) {
      warning_msg <- paste("Multiple matches for", role, "- columns:", paste(colnames(df)[match], collapse = ", "))
      warning(warning_msg)
      problems[[role]] <- warning_msg
      matched_cols[[role]] <- match[1]
    } else {
      matched_cols[[role]] <- match
    }
    cte <- cte + 1
  }
  col_roles <- unlist(matched_cols)
  if(anyDuplicated(col_roles, incomparables = NA))
    stop("Multiple column roles are assigned to the same column index.
         Please check for ambiguous or duplicate column matches.")

  # ---- 1. Extract matched columns safely ----
  snp_id_col     <- matched_cols[["snp_id"]]
  if(is.na(snp_id_col)) snp_id_col <- 1
  chrom_col <- matched_cols[["chrom"]]
  if(is.na(chrom_col)) chrom_col <- 2
  genome_pos_col <- matched_cols[["genome_pos"]]
  if(is.na(genome_pos_col)) genome_pos_col <- 3
  ref_col <- matched_cols[["ref"]]
  if(is.na(ref_col)) ref_col <- 4
  alt_col <- matched_cols[["alt"]]
  if(is.na(alt_col)) alt_col <- 5
  P1_col <- matched_cols[["P1"]]
  if(is.na(P1_col)) P1_col <- 6
  P2_col <- matched_cols[["P2"]]
  if(is.na(P2_col)) P2_col <- 7

  # ---- 2. Validate metadata columns ----
  if (!is.character(df[[snp_id_col]])) {
    warning_msg <- "SNP ID column must be of character type"
    warning(warning_msg)
    problems[["snp_id_type"]] <- warning_msg
  }
  if (anyDuplicated(df[[snp_id_col]])) {
    warning_msg <- "Duplicated SNP IDs found"
    warning(warning_msg)
    problems[["snp_id_duplicate"]] <- warning_msg
  }

  if(all(is.na(df[[chrom_col]]))){
    warning_msg <- "All values in the chromosome column are missing (NA)."
    warning(warning_msg)
    problems[["chrom_format"]] <- warning_msg
  } else if (!all(grepl("^(Chr_?\\d+|chr_?\\d+|Ch_?\\d+|ch_?\\d+|CHROM_?\\d+|CHR_?\\d+|\\d+|chr\\w+)$",
                      df[[chrom_col]]))) {
    warning_msg <- "Chromosome names must follow a numeric pattern\n   Usually Chr01, Chr02, etc"
    warning(warning_msg)
    problems[["chrom_format"]] <- warning_msg
  }

  # ---- 3. Check and optionally reorder by chromosome and position ----
  if(all(is.na(df[[genome_pos_col]]))){
    warning_msg <- "All values in the genome position column are missing (NA)."
    warning(warning_msg)
    problems[["genome_pos_type"]] <- warning_msg
  } else if (!is.numeric(df[[genome_pos_col]])) {
    warning_msg <- "Genome position column must be numeric"
    warning(warning_msg)
    problems[["genome_pos_type"]] <- warning_msg
  } else {
    unordered <- unlist(tapply(df[[genome_pos_col]], df[[chrom_col]], function(x) any(diff(x) < 0)))
    if (!reorder & any(unordered, na.rm = TRUE)) {
      warning_msg <- "Genome positions are not strictly increasing within some chromosomes"
      warning(warning_msg)
      problems[["genome_pos_order"]] <- warning_msg
    }
    if (reorder & any(unordered, na.rm = TRUE)) {
      warning_msg <- "Genome positions are not strictly increasing within some chromosomes\n   Data frame was reordered by chromosome and genome position."
      warning(warning_msg)
      problems[["genome_pos_order"]] <- warning_msg
      chr_order <- suppressWarnings(as.numeric(gsub("[^0-9]", "", df[[chrom_col]])))
      chr_order[is.na(chr_order)] <- Inf
      df <- df[order(chr_order, df[[genome_pos_col]]), , drop = FALSE]
    }
  }

  # ---- 4. Check ref and alt ----
  valid_nt <- c("A", "T", "C", "G", "-", " ", "")

  # Warn if any ref alleles are neither A/T/C/G nor NA
  invalid_ref <- !(df[[ref_col]] %in% valid_nt | is.na(df[[ref_col]]))
  if (any(invalid_ref, na.rm = TRUE)) {
    warning_msg <- "'ref' column contains invalid nucleotides (not A/T/C/G/-/' '/'' or NA)"
    warning(warning_msg)
    problems[["ref_invalid"]] <- warning_msg
  }

  # Same for alt
  invalid_alt <- !(df[[alt_col]] %in% valid_nt | is.na(df[[alt_col]]))
  if (any(invalid_alt, na.rm = TRUE)) {
    warning_msg <- "'alt' column contains invalid nucleotides (not A/T/C/G or NA)"
    warning(warning_msg)
    problems[["alt_invalid"]] <- warning_msg
  }

  # Only compare ref == alt for rows where both are valid nucleotides
  ref_valid <- df[[ref_col]] %in% valid_nt
  alt_valid <- df[[alt_col]] %in% valid_nt
  same_alleles <- df[[ref_col]] == df[[alt_col]] & ref_valid & alt_valid
  if (any(same_alleles, na.rm = TRUE)) {
    warning_msg <- "Some SNPs have identical ref and alt alleles (among valid nucleotides)"
    warning(warning_msg)
    problems[["ref_equals_alt"]] <- warning_msg
  }

  # Only perform consistency check for rows where at least one allele is not NA
  refalt_non_na <- !(is.na(df[[ref_col]]) & is.na(df[[alt_col]]))
  if (any(refalt_non_na)) {
    refalt_check <- df[refalt_non_na, ] %>%
      group_by(snp = .data[[colnames(df)[snp_id_col]]]) %>%
      summarise(
        n_ref = n_distinct(.data[[colnames(df)[ref_col]]], na.rm = TRUE),
        n_alt = n_distinct(.data[[colnames(df)[alt_col]]], na.rm = TRUE)
      ) %>%
      filter(n_ref > 1 | n_alt > 1)

    if (nrow(refalt_check) > 0) {
      warning_msg <- "Inconsistent ref/alt alleles found for some SNPs (excluding NAs)"
      warning(warning_msg)
      problems[["refalt_inconsistent"]] <- warning_msg
    }
  }

  # ---- 5. Dosage checks ----
  if (!is.numeric(df[[P1_col]]) || !is.numeric(df[[P2_col]])) {
    warning_msg <- paste0("P1 and P2 dosage columns (columns ", P1_col, " and ", P2_col, ") must be numeric")
    warning(warning_msg)
    problems[["dosage_type"]] <- warning_msg
  } else {
    if (any(df[[P1_col]] > ploidy.p1, na.rm = TRUE)) {
      warning_msg <- "Some P1 dosage values exceed ploidy.p1"
      warning(warning_msg)
      problems[["p1_exceeds_ploidy"]] <- warning_msg
    }
    if (any(df[[P2_col]] > ploidy.p2, na.rm = TRUE)) {
      warning_msg <- "Some P2 dosage values exceed ploidy.p2"
      warning(warning_msg)
      problems[["p2_exceeds_ploidy"]] <- warning_msg
    }
  }

  # ---- 6. Genotype checks ----
  metadata_cols <- c(snp_id_col, chrom_col, genome_pos_col, ref_col, alt_col)
  genotype_cols <- setdiff(seq_along(df), c(metadata_cols, 2, 3))
  if (any(!sapply(df[genotype_cols], is.numeric))) {
    warning_msg <- "Some genotype columns are not numeric"
    warning(warning_msg)
    problems[["genotype_type"]] <- warning_msg
  }

  max_ploidy <- max(ploidy.p1, ploidy.p2)
  if (any(unlist(df[genotype_cols]) > max_ploidy, na.rm = TRUE)) {
    warning_msg <- "Some genotype dosages exceed max ploidy"
    warning(warning_msg)
    problems[["genotype_exceeds_ploidy"]] <- warning_msg
  }

  # ---- 7. Missing data checks ----
  marker_na_rate <- rowMeans(is.na(df[genotype_cols]))
  if (any(marker_na_rate > max_missing_rate)) {
    warning_msg <- "Some SNPs have missing genotype rates > threshold"
    warning(warning_msg)
    problems[["na_rate"]] <- warning_msg
  }

  # ---- 8. Sample name checks ----
  sample_ids <- colnames(df)[genotype_cols]
  if (any(!grepl("^[A-Za-z0-9_.-]+$", sample_ids))) {
    warning_msg <- "Some sample IDs contain invalid characters"
    warning(warning_msg)
    problems[["sample_id_format"]] <- warning_msg
  }

  invisible(list(
    valid = TRUE,
    n_snps = nrow(df),
    n_samples = length(genotype_cols),
    ploidy.p1 = ploidy.p1,
    ploidy.p2 = ploidy.p2,
    max_na_rate_observed = max(marker_na_rate),
    problems = problems,
    df = df[,c(c(snp_id_col,
                 chrom_col,
                 genome_pos_col,
                 ref_col,
                 alt_col,
                 P1_col,
                 P2_col) , 8:ncol(df)), drop = FALSE]
  ))
}

# Function to check if an object is of class "mappoly2.data"
is.mappoly2.data <- function(x) {
  inherits(x, "mappoly2.data")
}

# Function to check if an object has been screened
has.mappoly2.screened <- function(x) {
  is.mappoly2.data(x) && length(x$screened.data) > 0
}

# Function to check if chromosome information is present
has.chromosome.info <- function(x) {
  is.mappoly2.data(x) &&
    !all(is.na(x$chrom)) && !is.null(x$chrom)
}

# Function to check if recombination frequency data is present
has.mappoly2.rf <- function(x) {
  has.mappoly2.screened(x) && inherits(x, "pairwise.rf")
}

# Function to check if genome position information is present
data.has.genome.info <- function(x) {
  is.mappoly2.data(x) && !all(is.na(x$genome.pos))
}

# Function to check if an object is a "mappoly2.sequence"
is.mappoly2.sequence <- function(x) {
  inherits(x, "mappoly2.sequence")
}

# Function to check if a sequence is phased
is.phased.sequence <- function(x) {
  is.mappoly2.sequence(x) &&
    !is.null(x$phases[[1]]$p1) &&
    !is.null(x$phases[[1]]$p2)
}

# Function to check if a sequence is phased based on specific parameters
is.phased.sequence <- function(x, lg, type, parent, phase.type) {
  assertthat::assert_that(length(lg) == 1)
  assertthat::assert_that(length(type) == 1)
  assertthat::assert_that(length(parent) == 1)
  assertthat::assert_that(length(phase.type) == 1)
  !is.null(x$maps[[lg]][[type]][[parent]][[phase.type]][[1]]$p1)
}

# Function to check if a sequence is MDS ordered
is.mds.ordered <- function(x, lg) {
  assertthat::assert_that(length(lg) == 1)
  !is.null(x$maps[[lg]][["mds"]]$order)
}

# Function to check if a sequence is mapped
is.mapped.sequence <- function(x, lg, type, parent) {
  assertthat::assert_that(length(lg) == 1)
  assertthat::assert_that(length(type) == 1)
  assertthat::assert_that(length(parent) == 1)
  !is.null(x$maps[[lg]][[type]][[parent]]$hmm.phase[[1]]$loglike)
}

# Function to check if a sequence is a haplotype sequence
is.haplotype.sequence <- function(x, lg, type, parent) {
  assertthat::assert_that(length(lg) == 1)
  assertthat::assert_that(length(type) == 1)
  assertthat::assert_that(length(parent) == 1)
  !is.null(x$maps[[lg]][[type]][[parent]]$hmm.phase[[1]]$haploprob)
}

# Function to check if ploidy levels are adequate (even)
has.adequate.ploidy <- function(x) {
  assertthat::assert_that(
    (x %% 2) == 0,
    msg = "At least one of the parents has an odd ploidy level."
  )
}

#' @title Check if SNPs are ordered by chromosome and position
#'
#' @description
#' Returns TRUE if chromosomes are in numeric order (e.g., Chr_1, Chr_2, ...)
#' and SNP positions are in ascending order within each chromosome.
#'
#' @param chrom A character vector of chromosome names.
#' @param pos A numeric vector of SNP positions.
#'
#' @return Logical. TRUE if both chromosome and position order is correct.
#'
#' @export
check_snp_order <- function(chrom, pos) {
  chrom_numeric <- suppressWarnings(as.numeric(gsub("[^0-9]", "", chrom)))
  chrom_ordered <- all(diff(chrom_numeric) >= 0)
  pos_ordered_by_chrom <- tapply(pos, chrom, function(p) all(diff(p) >= 0))
  all_pos_ordered <- all(pos_ordered_by_chrom)
  list(ch_ord = chrom_ordered, snp_ord = all_pos_ordered)
}

