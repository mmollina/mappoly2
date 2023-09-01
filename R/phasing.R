#' Phasing based on pairwise recombination fraction estimation
#'
#' @param void internal function
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
pairwise_phasing_per_group <- function(x,
                                       gr,
                                       type = c("mds", "genome"),
                                       thresh.LOD.ph = 3,
                                       thresh.LOD.rf = 3,
                                       thresh.rf = 0.5,
                                       max.conf.btnk.p1 = 1,
                                       max.conf.btnk.p2 = max.conf.btnk.p1,
                                       verbose = TRUE){
  assert_that(inherits(x, "ws"))
  type<- match.arg(type)
  if(type == "mds")
    mrk.id <- as.character(x$working.sequences[[gr]]$order$mds$info$locimap$locus)
  if(type == "genome")
    mrk.id <- as.character(x$working.sequences[[gr]]$order$genome$info$mrk.names) #FIXME
  m<-mappoly2:::filter_rf_matrix(x,
                                 type = "sh",
                                 thresh.LOD.ph,
                                 thresh.LOD.rf,
                                 thresh.rf)
  assert_that(all(rownames(mrk.id)%in%x$pairwise$rec.mat), msg = "some markers in the sequence are not present in the pairwise rf.")
  {
    if(verbose)
      cat("Phasing parent", x$data$name.p1, "\n")
    Ph.p1 <- mappoly2:::twopt_phasing_cpp(mrk_id = mrk.id,
                                          ploidy = x$data$ploidy.p1,
                                          dose_vec = x$data$dosage.p1[mrk.id],
                                          S = m$Sh.p1[mrk.id, mrk.id],
                                          max_conf_number = max.conf.btnk.p1,
                                          verbose = verbose)
    for(i in 1:length(Ph.p1$phase_configs))
      rownames(Ph.p1$phase_configs[[i]]) <- Ph.p1$marker_names
    if(verbose)
      cat("Phasing parent", x$data$name.p2, "\n")
    Ph.p2 <- mappoly2:::twopt_phasing_cpp(mrk_id = mrk.id,
                                          ploidy = x$data$ploidy.p2,
                                          dose_vec = x$data$dosage.p2[mrk.id],
                                          S = m$Sh.p2[mrk.id, mrk.id],
                                          max_conf_number = max.conf.btnk.p2,
                                          verbose = verbose)
    for(i in 1:length(Ph.p2$phase_configs))
      rownames(Ph.p2$phase_configs[[i]]) <- Ph.p2$marker_names
  }
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
                        p2 = Ph.p2$phase_configs[[j]][mrks,],
                        loglike = NULL,
                        rf = NULL,
                        error = NULL,
                        haploprob = NULL)
      cte <- cte + 1
    }
  }
  if(type == "mds")
    x$working.sequences[[gr]]$order$mds$phase = unique(Ph)
  if(type == "genome")
    x$working.sequences[[gr]]$order$genome$phase = unique(Ph)
  return(x)
}




