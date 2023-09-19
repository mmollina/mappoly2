#' Conditional probabilities computation using Hidden Markov Models
#'
#' @param void internal function
#' @keywords internal
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
calc_haplotypes_per_group <- function(x,
                                      gr,
                                      type = c("both", "mds", "genome"))
{
  type <- match.arg(type)
  id.map <- sapply(x$working.sequences$Lg_1$order,
                   function(y) !is.null(y$phase[[1]]$loglike))
  if(all(!id.map)) stop("the map is not estimated")

  for(i in which(id.map)){

  }




  mrk.id <- rownames(x$working.sequencesphases[[1]]$p1)
  g <- x$data$geno.dose[mrk.id, ]
  g[is.na(g)] <- -1
  phase.conf <- match.arg(phase.conf)

  if(phase.conf == "all")
    phase.conf <- 1:length(x$phases)

  if(phase.conf == "best")
    phase.conf <- 1

  output.seq <- x

  cat("Computing haplotype probabilities\n")
  cat("   Number of phase configurations: ", length(phase.conf), "\n")
  if (detect_info_par(x) == "both" || compute.both.parents){
    for(i in 1: length(phase.conf)){
      cat("   Conf.", i,":")
      pedigree <- matrix(rep(c(1,
                               2,
                               x$data$ploidy.p1,
                               x$data$ploidy.p2, 1),
                             x$data$n.ind),
                         nrow = x$data$n.ind,
                         byrow = TRUE)
      output.seq$phases[[phase.conf[i]]]$haploprob <- calc_haploprob_biallelic(PH = list(x$phases[[i]]$p1,
                                                                                        x$phases[[i]]$p2),
                                                                              G = g,
                                                                              pedigree = pedigree,
                                                                              rf = x$phases[[phase.conf[i]]]$rf,
                                                                              err = x$phases[[phase.conf[i]]]$error)
      cat("\n")
    }
    cat("Done with haplotype probabilities\n")
    return(output.seq)
  } else if(detect_info_par(x) == "p1"){
    id <- which(x$data$ploidy.p2 == x$data$dosage.p2[mrk.id])
    g[id, ] <- g[id, ] - x$data$ploidy.p2/2
    for(i in phase.conf){
      cat("   Conf.", i,":")
      output.seq$phases[[phase.conf[i]]]$haploprob <- calc_haploprob_biallelic_single(PH = x$phases[[i]]$p1,
                                                                                     G = g,
                                                                                     rf = x$phases[[phase.conf[i]]]$rf,
                                                                                     err = x$phases[[phase.conf[i]]]$error)
      cat("\n")
    }
    cat("Done with haplotype probabilities\n")
    return(output.seq)
  } else if(detect_info_par(x) == "p2") {
    id <- which(x$data$ploidy.p1 == x$data$dosage.p1[mrk.id])
    g[id, ] <- g[id, ] - x$data$ploidy.p1/2
    for(i in phase.conf){
      cat("   Conf.", i,":")
      output.seq$phases[[phase.conf[i]]]$haploprob <- calc_haploprob_biallelic_single(PH = x$phases[[i]]$p2,
                                                                                   G = g,
                                                                                   rf = x$phases[[phase.conf[i]]]$rf,
                                                                                   err = x$phases[[phase.conf[i]]]$error)
    }
    cat("Done with haplotype probabilities\n")
    return(output.seq)
  } else {
    stop("it should not get here")
  }
}
