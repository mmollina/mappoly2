#' Multipoint analysis using Hidden Markov Models
#'
#' @param void internal function
#' @keywords internal
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
hmm_map_reconstruction <- function(input.seq,
                                   phase.conf = "all",
                                   rf = NULL,
                                   verbose = FALSE,
                                   tol = 10e-4)
{
  assert_that(is.mappoly2.sequence(input.seq))
  mrk.id <- rownames(input.seq$phases[[1]]$p1)
  g <- input.seq$data$geno.dose[mrk.id, ]
  ret_H0 <- TRUE
  if(is.null(rf)){
    rf <- rep(0.01, nrow(g) - 1)
    ret_H0 <- FALSE
  }
  assert_that(length(rf) == nrow(g) - 1)
  if(phase.conf == "all")
    phase.conf <- 1:length(input.seq$phases)

  output.seq <- input.seq

  assert_that(all(phase.conf%in%1:length(input.seq$phases)),
              msg = "invalid phases specified in 'phase.conf'")
  cat("Multi-locus map estimation\n")
  cat("   Number of phase configurations: ", length(phase.conf), "\n")
  if (detect_info_par(input.seq) == "both"){
    for(i in phase.conf){
      cat("   Conf.", i,":")
      pedigree <- matrix(rep(c(1,
                               2,
                               input.seq$data$ploidy.p1,
                               input.seq$data$ploidy.p2, 1),
                             input.seq$data$n.ind),
                         nrow = input.seq$data$n.ind,
                         byrow = TRUE)
      w <- est_hmm_map_biallelic(PH = list(input.seq$phases[[i]]$p1,
                                           input.seq$phases[[i]]$p2),
                                 G = g,
                                 pedigree = pedigree,
                                 rf = rf,
                                 verbose = verbose,
                                 detailed_verbose = FALSE,
                                 tol = tol,
                                 ret_H0 = ret_H0)
      output.seq$phases[[i]]$loglike <- w[[1]]
      output.seq$phases[[i]]$rf <- w[[2]]
    }
    cat("Done with map estimation ~~~~~~ \n")
    return(output.seq)
  } else if(detect_info_par(input.seq) == "p1"){
    id <- which(input.seq$data$ploidy.p2 == input.seq$data$dosage.p2[mrk.id])
    g[id, ] <- g[id, ] - input.seq$data$ploidy.p2/2
    for(i in phase.conf){
      cat("   Conf.", i,":")
      w <- est_hmm_map_biallelic_single(PH = input.seq$phases[[i]]$p1,
                                        G = g,
                                        rf = rf,
                                        verbose = verbose,
                                        detailed_verbose = FALSE,
                                        tol = tol,
                                        ret_H0 = ret_H0)
      output.seq$phases[[i]]$loglike <- w[[1]]
      output.seq$phases[[i]]$rf <- w[[2]]
    }
    cat("Done with map estimation\n")
    return(output.seq)
  } else if(detect_info_par(input.seq) == "p2") {
    id <- which(input.seq$data$ploidy.p1 == input.seq$data$dosage.p1[mrk.id])
    g[id, ] <- g[id, ] - input.seq$data$ploidy.p1/2
    for(i in phase.conf){
      cat("   Conf.", i,":")
      w <- est_hmm_map_biallelic_single(PH = input.seq$phases[[i]]$p2,
                                        G = g,
                                        rf = rf,
                                        verbose = verbose,
                                        detailed_verbose = FALSE,
                                        tol = tol,
                                        ret_H0 = ret_H0)
      output.seq$phases[[i]]$loglike <- w[[1]]
      output.seq$phases[[i]]$rf <- w[[2]]
    }
    cat("Done with map estimation\n")
    return(output.seq)
  } else {
    stop("it should not get here")
  }
}
