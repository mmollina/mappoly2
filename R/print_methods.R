#' @export
print.mappoly2.data <- function(x, type = c("screened", "raw"), detailed = FALSE,  ...) {
  #### Raw ####
  msg("Data summary", col = "blue")
  txt <- list(
    paste0("    Ploidy level of ", x$name.p1, ":"),
    paste0("    Ploidy level of ", x$name.p2, ":"),
    "    No. individuals:",
    "    No. markers:",
    "    Percentage of missing:",
    "    Chromosome info:",
    "    Genome position:",
    "    Recombination farction:",
    "    Marker scope:")
  n <- sapply(txt, nchar)
  for (i in 1:length(txt)) {
    txt[[i]] <- paste(txt[[i]], paste0(rep(" ", max(n) - n[i]), collapse = ""))
  }
  id <- is.na(x$geno.dose)
  cat(" ", txt[[1]], " ", x$ploidy.p1, sep ="")
  cat("\n", txt[[2]], x$ploidy.p2)
  cat("\n", txt[[3]], x$n.ind)
  cat("\n", txt[[4]], length(x$mrk.names))
  cat("\n ", txt[[5]], " ",   round(100*sum(id)/length(id),1), "%" ,sep = "")
  chrom.flag <- FALSE
  if (all(is.null(x$chrom)) || all(is.na(x$chrom)))
    cat("\n", txt[[6]], "unavailable")
  else{
    cat("\n", txt[[6]], "available")
    chrom.flag <- TRUE
  }
  if(any(is.numeric(x$genome.pos)))
    cat("\n",  txt[[7]], "available")
  else
    cat("\n",  txt[[7]], "unavailable")

  #### RF ####
    if(has.mappoly2.rf(x)){
      cat("\n",  txt[[8]], "available")
      cat("\n ",  txt[[9]], " ",
          x$pairwise.rf$mrk.scope,
          " (",
          paste(unique(x$chrom[colnames(x$pairwise.rf$rec.mat)]), collapse = ", "),
          ")\n\n", sep = "")
    } else
    cat("\n",  txt[[8]], "unavailable\n\n")
  #### Screened ####
  if(mappoly2:::has.mappoly2.screened(x)){
    msg("Filtering information", col = "blue")
    txt <- list(
      "    Thresholds",
      "       missing mrk: ",
      "       missing ind: ",
      "       chi-square pval: ",
      "       read depth: ",
      "       non full-sib: ",
      "    Screened mrk: ",
      "    Screened ind: ")
    n <- sapply(txt, nchar)
    for (i in 1:length(txt)) {
      txt[[i]] <- paste(txt[[i]], paste0(rep(" ", max(n) - n[i]), collapse = ""))
    }
    id <- is.na(x$geno.dose)
    cat(txt[[1]])
    cat("\n", txt[[2]], x$screened.data$thresholds$miss.mrk)
    cat("\n", txt[[3]], x$screened.data$thresholds$miss.ind)
    cat("\n", txt[[4]], format(x$screened.data$thresholds$chisq.pval, digits = 3))
    cat("\n", txt[[5]], x$screened.data$thresholds$read.depth)
    cat("\n", txt[[6]], ifelse(all(is.na(x$QAQC.values$individuals$full.sib)),
                               "-", sum(!x$QAQC.values$individuals$full.sib)))
    cat("\n", txt[[7]], length(x$screened.data$mrk.names))
    cat("\n", txt[[8]], length(x$screened.data$ind.names), "\n\n")
  }
  #### Detailed ####
  w <- table(x$chrom, useNA = "always")
  w <- w[order(mappoly2:::embedded_to_numeric(names(w)))]
  names(w)[is.na(names(w))] <- "NoCrh"
  if(detailed){
    if(chrom.flag){
      msg("No. markers per chromosome", col = "blue")
      print(data.frame(chrom = paste0("       ", names(w)), No.mrk = as.numeric(w)), row.names = FALSE)
    }
    cat("\n")
    msg("No. of markers per dosage in both parents", col = "blue")
    freq <- table(paste(x$dosage.p1,
                        x$dosage.p2, sep = "-"))
    d.temp <- matrix(unlist(strsplit(names(freq), "-")), ncol = 2, byrow = TRUE)
    d.temp <- data.frame(paste0("    ", d.temp[, 1]),
                         d.temp[, 2],
                         as.numeric(freq))
    colnames(d.temp) <- c(x$name.p1, x$name.p2, "freq")
    print(d.temp, row.names = FALSE)
  }
}
