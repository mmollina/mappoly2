#' Phasing Based on Pairwise Recombination Fraction Estimation
#'
#' This function performs phasing based on pairwise recombination fraction estimation for genetic maps in a mappoly2 sequence object. It is designed to facilitate the phasing process using detailed genetic data.
#'
#' @param x An object representing genetic mapping data, typically a mappoly2 sequence object.
#' @param lg Optional vector specifying the linkage groups to be processed.
#'           If NULL, all linkage groups in `x` are considered.
#' @param type The type of genetic maps to be processed, options include "mds", "genome", or "custom".
#' @param parent Specifies which parent's data to use in the phasing process.
#'               Options are "p1p2" (both parents), "p1" (first parent), or "p2" (second parent).
#' @param thresh.LOD.ph Threshold for the LOD (Logarithm of the Odds) score for phasing.
#' @param thresh.LOD.rf Threshold for the LOD score for recombination fractions.
#' @param thresh.rf Threshold for recombination fraction.
#' @param max.search.expansion.p1 The maximum number of search expansions for parent 1.
#' @param max.search.expansion.p2 The maximum number of search expansions for parent 2.
#'                                Defaults to the same as `max.search.expansion.p1`.
#' @param verbose Logical; if TRUE, progress messages will be printed.
#'
#' @return Returns the updated mappoly2 sequence object with the phased genetic maps.
#'
#' @details The function iterates over the specified linkage groups and performs phasing
#'          based on pairwise recombination fraction estimation. It supports phasing for
#'          individual parents as well as both parents together. The function handles the
#'          alignment and integration of phasing data across different maps.
#'
#' @importFrom assertthat assert_that
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
pairwise_phasing <- function(x,
                             lg = NULL,
                             type = c("mds", "genome"),
                             parent = c("p1p2","p1","p2"),
                             thresh.LOD.ph = 3,
                             thresh.LOD.rf = 3,
                             thresh.rf = 0.5,
                             max.search.expansion.p1 = 10,
                             max.search.expansion.p2 = max.search.expansion.p1,
                             verbose = TRUE){
  y <- parse_lg_and_type(x,lg,type)
  assert_that(has.mappoly2.screened(x$data))
  parent <- match.arg(parent)
  mrk.id <- get_markers_from_ordered_sequence(x, y$lg, y$type, parent)
  for(i in 1:length(mrk.id)){
    if(verbose) cat("  -->", y$lg[i], "\n")
    x$maps[[y$lg[i]]][[y$type]][[parent]]$rf.phase <- pairwise_phasing_one(x,
                                                                           mrk.id[[i]],
                                                                           thresh.LOD.ph,
                                                                           thresh.LOD.rf,
                                                                           thresh.rf,
                                                                           max.search.expansion.p1,
                                                                           max.search.expansion.p2,
                                                                           verbose)
    if(parent == "p1p2"){
      ## attributing phase to p1 informative markers
      temp <- x$maps[[y$lg[i]]][[y$type]][[parent]]$rf.phase
      idx <- get_info_markers(rownames(temp[[1]]$p1), x, "p1")
      for(j in 1:length(temp)){
        temp[[j]]$p1 <- temp[[j]]$p1[idx,]
        temp[[j]]$p2 <- temp[[j]]$p2[idx,]
      }
      x$maps[[y$lg[i]]][[y$type]][["p1"]]$rf.phase <- temp
      ## attributing phase to p1 informative markers
      temp <- x$maps[[y$lg[i]]][[y$type]][[parent]]$rf.phase
      idx <- get_info_markers(rownames(temp[[1]]$p1), x, "p2")
      for(j in 1:length(temp)){
        temp[[j]]$p1 <- temp[[j]]$p1[idx,]
        temp[[j]]$p2 <- temp[[j]]$p2[idx,]
      }
      x$maps[[y$lg[i]]][[y$type]][["p2"]]$rf.phase <- temp
    }
    if(verbose) cat("\n")
  }
  return(x)
}

pairwise_phasing_one <- function(x,
                                 mrk.id,
                                 thresh.LOD.ph = 3,
                                 thresh.LOD.rf = 3,
                                 thresh.rf = 0.5,
                                 max.search.expansion.p1,
                                 max.search.expansion.p2 = max.search.expansion.p1,
                                 verbose = TRUE){
  assert_that(!is.null(mrk.id), msg = "Provide an ordered sequence.")
  assert_that(has.mappoly2.screened(x$data))
  assert_that(all(mrk.id %in% x$data$screened.data$mrk.names))
  M <- filter_rf_matrix(x$data,
                        type = "sh",
                        thresh.LOD.ph,
                        thresh.LOD.rf,
                        thresh.rf,
                        mrk.names = mrk.id)
  if(verbose)
    cat("Phasing parent", x$data$name.p1, "\n")
  Ph.p1 <- twopt_phasing_cpp(mrk_id = mrk.id,
                             ploidy = x$data$ploidy.p1,
                             dose_vec = x$data$dosage.p1[mrk.id],
                             S = M$Sh.p1,
                             max_conf_number = max.search.expansion.p1,
                             verbose = verbose)
  for(i in 1:length(Ph.p1$phase_configs))
    rownames(Ph.p1$phase_configs[[i]]) <- Ph.p1$marker_names
  if(verbose)
    cat("Phasing parent", x$data$name.p2, "\n")
  Ph.p2 <- twopt_phasing_cpp(mrk_id = mrk.id,
                             ploidy = x$data$ploidy.p2,
                             dose_vec = x$data$dosage.p2[mrk.id],
                             S = M$Sh.p2,
                             max_conf_number = max.search.expansion.p2,
                             verbose = verbose)
  for(i in 1:length(Ph.p2$phase_configs))
    rownames(Ph.p2$phase_configs[[i]]) <- Ph.p2$marker_names
  mrks <- intersect(Ph.p1$marker_names, Ph.p2$marker_names)
  n1 <- length(mrk.id)
  n2 <- length(mrks)
  if(verbose){
    cat(n2, " phased markers out of ", n1, ": (",  round(100*n2/n1,1), "%)",sep = "")
  }
  cte <- 1
  Ph <- vector("list", length(Ph.p1$phase_configs) * length(Ph.p2$phase_configs))
  for(i in 1:length(Ph.p1$phase_configs)){
    for(j in 1:length(Ph.p2$phase_configs)){
      Ph[[cte]] <- list(p1 = Ph.p1$phase_configs[[i]][mrks,],
                        p2 = Ph.p2$phase_configs[[j]][mrks,])
      cte <- cte + 1
    }
  }
  return(unique(Ph))
}
