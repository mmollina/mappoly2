#' Export biparental genetic map to CSV
#'
#' This function exports the mapped marker information from a biparental population
#' stored in a `mappoly2.sequence` object to a CSV file. It includes linkage group,
#' physical position, allele information, dosages, and parental haplotypes.
#'
#' @param x A `mappoly2.sequence` object containing phased maps.
#' @param type Type of map to export. Must be either `"mds"` or `"genome"`. Default is `"mds"`.
#' @param parent Parent phase to use: `"p1p2"`, `"p1"`, or `"p2"`. Default is `"p1p2"`.
#' @param file Output file name. Default is `"map_output.csv"`.
#'
#' @return Invisibly returns a data frame containing the map that was written to the CSV file.
#' @export
#'
#' @examples
#' \dontrun{
#'   export_biparental_map_to_csv(my_seq_object, type = "genome", parent = "p1p2", file = "my_map.csv")
#' }
export_biparental_map_to_csv <- function(x,
                                         type = c("mds", "genome"),
                                         parent = c("p1p2", "p1", "p2"),
                                         file = "map_output.csv") {
  assertthat::assert_that(is.mappoly2.sequence(x))

  type <- match.arg(type)
  parent <- match.arg(parent)

  # Check if map is available for each linkage group
  has_map <- vapply(names(x$maps), function(i) {
    mappoly2:::is.mapped.sequence(x, i, type, parent)
  }, logical(1))
  names(has_map) <- names(x$maps)

  if (!all(has_map)) {
    stop("Maps are not available for all linkage groups for:\n--> type: ", type, "\n--> parent: ", parent)
  }

  results <- list()
  for (lg in names(x$maps)[has_map]) {
    map <- x$maps[[lg]][[type]]
    phase_info <- map[[parent]]$hmm.phase[[1]]
    ph_p1 <- phase_info$p1
    ph_p2 <- phase_info$p2
    marker_names <- rownames(ph_p1)

    # Replace 0/1 by reference/alternate alleles
    for (j in seq_along(marker_names)) {
      ref <- x$data$ref[marker_names[j]]
      alt <- x$data$alt[marker_names[j]]
      ph_p1[j, ] <- ifelse(ph_p1[j, ] == 0, ref, alt)
      ph_p2[j, ] <- ifelse(ph_p2[j, ] == 0, ref, alt)
    }

    ph_p1 <- as.data.frame(ph_p1)
    ph_p2 <- as.data.frame(ph_p2)
    colnames(ph_p1) <- paste0("P1_", seq_len(x$data$ploidy.p1))
    colnames(ph_p2) <- paste0("P2_", seq_len(x$data$ploidy.p2))

    df <- dplyr::tibble(
      `Marker Name` = marker_names,
      `LG` = rep(lg, length(marker_names)),
      `Ref Chrom` = x$data$chrom[marker_names],
      `Ref Position` = x$data$genome.pos[marker_names],
      `Ref Allele` = x$data$ref[marker_names],
      `Alt Allele` = x$data$alt[marker_names],
      `Map Position` = round(cumsum(imf_h(c(0, phase_info$rf))), 2),
      `Dose in P1` = x$data$dosage.p1[marker_names],
      `Dose in P2` = x$data$dosage.p2[marker_names]
    ) |>
      dplyr::bind_cols(ph_p1, ph_p2)

    results[[lg]] <- df
  }

  final_df <- do.call(rbind, results)
  write.csv(final_df, file = file, row.names = FALSE)
  invisible(final_df)
}

#' Export a `mappoly2.data` Object to a CSV File
#'
#' This function exports a \code{mappoly2.data} object to a CSV file containing
#' SNP metadata and dosage information for both parents and all offspring.
#'
#' @param x An object of class \code{mappoly2.data}.
#' @param file A character string naming the output CSV file.
#' @param name.p1 Optional name to replace the default label \code{P1} (parent 1).
#' @param name.p2 Optional name to replace the default label \code{P2} (parent 2).
#'
#' @return Invisibly returns a \code{data.frame} containing the exported data.
#'
#' @export
export_data_to_csv <- function(x, file = "file.csv",
                               name.p1 = NULL,
                               name.p2 = NULL) {
  assertthat::assert_that(inherits(x, "mappoly2.data"))

  # Use custom parent names if provided
  p1_name <- if (!is.null(name.p1)) name.p1 else "P1"
  p2_name <- if (!is.null(name.p2)) name.p2 else "P2"

  # Build output data frame
  out <- tibble::tibble(
    snp_id     = x$mrk.names,
    Chr        = x$chrom,
    genome_pos = x$genome.pos,
    ref        = x$ref,
    alt        = x$alt
  )
  out[[p1_name]] <- x$dosage.p1
  out[[p2_name]] <- x$dosage.p2

  # Add progeny dosages
  geno <- as.data.frame(x$geno.dose)
  out <- cbind(out, geno)

  # Write to CSV
  write.csv(out, file = file, row.names = FALSE)

  invisible(out)
}

