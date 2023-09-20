#' Conditional probabilities computation using Hidden Markov Models
#'
#' @param void internal function
#' @keywords internal
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
calc_haplotypes_per_group <- function(x, gr, type = c("both", "mds", "genome"))
{
  type <- match.arg(type)
  assert_that(inherits(x, "ws"))
  idx <- sapply(x$working.sequences[[gr]]$order, function(x) !is.null(x$hmm.map$loglike))
  if(type!="both")
    idx <- idx[type]
  assert_that(any(idx))
  cat("Computing haplotype probabilities\n")
  cat("  Sequence: ", gr, "\n")
  for(i in which(idx)){
    cat("   '-->", names(idx)[i])
    mrk.id <- rownames(x$working.sequences[[gr]]$order[[i]]$hmm.map$p1)
    ind.id <- x$working.sequences[[gr]]$ind.names
    n.ind <- length(ind.id)
    g <- x$data$geno.dose[mrk.id, ind.id]
    g[is.na(g)] <- -1
    if (detect_info_par(x, gr) == "both"){
      pedigree <- matrix(rep(c(1,
                               2,
                               x$data$ploidy.p1,
                               x$data$ploidy.p2, 1),
                             n.ind),
                         nrow = n.ind,
                         byrow = TRUE)
      x$working.sequences[[gr]]$order[[i]]$hmm.map$haploprob <- calc_haploprob_biallelic(PH = list(x$working.sequences[[gr]]$order[[i]]$hmm.map$p1,
                                                                                                 x$working.sequences[[gr]]$order[[i]]$hmm.map$p2),
                                                                                       G = g,
                                                                                       pedigree = pedigree,
                                                                                       rf = x$working.sequences[[gr]]$order[[i]]$hmm.map$rf,
                                                                                       err = x$working.sequences[[gr]]$order[[i]]$hmm.map$error)
    } else if(detect_info_par(x) == "p1"){
      id <- which(x$data$ploidy.p2 == x$data$dosage.p2[mrk.id])
      g[id, ] <- g[id, ] - x$data$ploidy.p2/2
      x$working.sequences[[gr]]$order[[i]]$hmm.map$haploprob <- calc_haploprob_biallelic_single(PH = x$working.sequences[[gr]]$order[[i]]$hmm.map$p1,
                                                                                              G = g,
                                                                                              rf = x$working.sequences[[gr]]$order[[i]]$hmm.map$rf,
                                                                                              err = x$working.sequences[[gr]]$order[[i]]$hmm.map$error)
    } else if(detect_info_par(x) == "p2"){
      id <- which(x$data$ploidy.p1 == x$data$dosage.p1[mrk.id])
      g[id, ] <- g[id, ] - x$data$ploidy.p1/2
      x$working.sequences[[gr]]$order[[i]]$hmm.map$haploprob <- calc_haploprob_biallelic_single(PH = x$working.sequences[[gr]]$order[[i]]$hmm.map$p2,
                                                                                              G = g,
                                                                                              rf = x$working.sequences[[gr]]$order[[i]]$hmm.map$rf,
                                                                                              err = x$working.sequences[[gr]]$order[[i]]$hmm.map$error)
    }
    cat("\n")
  }
  return(x)
}
