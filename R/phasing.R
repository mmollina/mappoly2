#' Phasing based on pairwise recombination fraction estimation
#'
#' @param void internal function
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
pairwise_phasing <- function(x,
                             lg = NULL,
                             type = c("mds", "genome"),
                             parent = c("p1p2","p1","p2"),
                             thresh.LOD.ph = 5,
                             thresh.LOD.rf = 5,
                             thresh.rf = 0.5,
                             max.search.expansion.p1 = 10,
                             max.search.expansion.p2 = max.search.expansion.p1,
                             verbose = TRUE){
  y <- mappoly2:::parse_lg_and_type(x,lg,type)
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
  assert_that(all(mrk.id %in% dat$screened.data$mrk.names))
  M <- mappoly2:::filter_rf_matrix(x$data,
                                   type = "sh",
                                   thresh.LOD.ph,
                                   thresh.LOD.rf,
                                   thresh.rf,
                                   mrk.names = mrk.id)
  if(verbose)
    cat("Phasing parent", x$data$name.p1, "\n")
  Ph.p1 <- mappoly2:::twopt_phasing_cpp(mrk_id = mrk.id,
                                        ploidy = x$data$ploidy.p1,
                                        dose_vec = x$data$dosage.p1[mrk.id],
                                        S = M$Sh.p1,
                                        max_conf_number = max.search.expansion.p1,
                                        verbose = verbose)
  for(i in 1:length(Ph.p1$phase_configs))
    rownames(Ph.p1$phase_configs[[i]]) <- Ph.p1$marker_names
  if(verbose)
    cat("Phasing parent", x$data$name.p2, "\n")
  Ph.p2 <- mappoly2:::twopt_phasing_cpp(mrk_id = mrk.id,
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
