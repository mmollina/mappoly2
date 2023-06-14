#' Conditional probabilities computation using Hidden Markov Models
#'
#' @param void internal function
#' @keywords internal
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
calc_haplotypes <- function(input.seq,
                            phase.conf = c("best","all"),
                            verbose = FALSE,
                            tol = 10e-4)
{
  assert_that(is.mappoly2.sequence(input.seq))
  mrk.id <- rownames(input.seq$phases[[1]]$p1)
  g <- input.seq$data$geno.dose[mrk.id, ]
  g[is.na(g)] <- -1
  phase.conf <- match.arg(phase.conf)

  if(phase.conf == "all")
    phase.conf <- 1:length(input.seq$phases)

  if(phase.conf == "best")
    phase.conf <- 1

  output.seq <- input.seq

  cat("Computing haplotype probabilities\n")
  cat("   Number of phase configurations: ", length(phase.conf), "\n")
  if (detect_info_par(input.seq) == "both"){
    for(i in 1: length(phase.conf)){
      cat("   Conf.", i,":")
      pedigree <- matrix(rep(c(1,
                               2,
                               input.seq$data$ploidy.p1,
                               input.seq$data$ploidy.p2, 1),
                             input.seq$data$n.ind),
                         nrow = input.seq$data$n.ind,
                         byrow = TRUE)
      output.seq$phases[[phase.conf[i]]]$haploprob <- calc_haploprob_biallelic(PH = list(input.seq$phases[[i]]$p1,
                                                                                        input.seq$phases[[i]]$p2),
                                                                              G = g,
                                                                              pedigree = pedigree,
                                                                              rf = input.seq$phases[[phase.conf[i]]]$rf,
                                                                              err = input.seq$phases[[phase.conf[i]]]$error)
      cat("\n")
    }
    cat("Done with haplotype probabilities\n")
    return(output.seq)
  } else if(detect_info_par(input.seq) == "p1"){
    id <- which(input.seq$data$ploidy.p2 == input.seq$data$dosage.p2[mrk.id])
    g[id, ] <- g[id, ] - input.seq$data$ploidy.p2/2
    for(i in phase.conf){
      cat("   Conf.", i,":")
      output.seq$phases[[phase.conf[i]]]$haploprob <- calc_haploprob_biallelic_single(PH = input.seq$phases[[i]]$p1,
                                                                                     G = g,
                                                                                     rf = input.seq$phases[[phase.conf[i]]]$rf,
                                                                                     err = input.seq$phases[[phase.conf[i]]]$error)
      cat("\n")
    }
    cat("Done with haplotype probabilities\n")
    return(output.seq)
  } else if(detect_info_par(input.seq) == "p2") {
    id <- which(input.seq$data$ploidy.p1 == input.seq$data$dosage.p1[mrk.id])
    g[id, ] <- g[id, ] - input.seq$data$ploidy.p1/2
    for(i in phase.conf){
      cat("   Conf.", i,":")
      output.seq$phases[[phase.conf[i]]]$haploprob <- calc_haploprob_biallelic_single(PH = input.seq$phases[[i]]$p2,
                                                                                   G = g,
                                                                                   rf = input.seq$phases[[phase.conf[i]]]$rf,
                                                                                   err = input.seq$phases[[phase.conf[i]]]$error)
    }
    cat("Done with haplotype probabilities\n")
    return(output.seq)
  } else {
    stop("it should not get here")
  }
}
