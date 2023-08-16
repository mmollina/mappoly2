#' @export
print.mappoly2 <- function(x, detailed = FALSE,  ...) {
  txt <- list(
    paste0("    Ploidy level of ", x$data$name.p1, ":"),
    paste0("    Ploidy level of ", x$data$name.p2, ":"),
    paste0("    No. individuals:"),
    paste0("    No. markers:"),
    paste0("    Percentage of missing:"),
    paste0("    Pairwise:"),
    paste0("    Linkage Groups:"),
    paste0("    Phases:"),
    paste0("      '--> Number of configurations:"),
    paste0("      '--> Percentage phased:"))
  n <- sapply(txt, nchar)
  for (i in 1:length(txt)) {
    txt[[i]] <- paste(txt[[i]], paste0(rep(" ", max(n) - n[i]), collapse = ""))
  }
  id <- is.na(x$data$geno.dose)
  cat("\n", txt[[1]], x$data$ploidy.p1)
  cat("\n", txt[[2]], x$data$ploidy.p2)
  cat("\n", txt[[3]], x$data$n.ind)
  cat("\n", txt[[4]], length(x$data$mrk.names))
  cat("\n ", txt[[5]], " ",   round(100*sum(id)/length(id),1), "%", sep = "")
  cat("\n", txt[[6]], ifelse(is.null(x$pairwise), "Not estimated", "Estimated"))
  cat("\n", txt[[7]])
  if(is.null(x$linkage.groups))
    cat(" Not allocated")
  else
    print_group(x$linkage.groups)
  cat("\n", txt[[8]])
  if(is.null(x$phases)){
    cat("\n ", txt[[9]], " 0", sep = "")
    cat("\n ", txt[[10]], " 0%", sep = "")
  } else {
    cat("\n ", txt[[9]], " ", length(x$phases), sep = "")
    cat("\n ", txt[[10]], " ", nrow(x$phases[[1]]$p1), " (",   round(100*nrow(x$phases[[1]]$p1)/length(x$mrk.names),1), "%)", sep = "")
  }
  w <- table(x$data$chrom, useNA = "always")
  w <- w[order(as.integer(gsub("[^0-9]", "", names(w))))]
  names(w)[is.na(names(w))] <- "NoCrh"
  if (all(is.null(x$data$chrom)) || all(is.na(x$data$chrom)))
    cat("\n     No. markers per chromosome: not available")
  else if(detailed){
    cat("\n     ----------\n     No. markers per chromosome:\n")
    print(data.frame(chrom = paste0("       ", names(w)), No.mrk = as.numeric(w)), row.names = FALSE)
  } else {
    cat("\n")
  }
    cat("     ----------\n     No. of markers per dosage in both parents:\n")
    freq <- table(paste(x$data$dosage.p1,
                        x$data$dosage.p2, sep = "-"))
    d.temp <- matrix(unlist(strsplit(names(freq), "-")), ncol = 2, byrow = TRUE)
    d.temp <- data.frame(paste0("    ", d.temp[, 1]),
                         d.temp[, 2],
                         as.numeric(freq))
    colnames(d.temp) <- c(x$data$name.p1, x$data$name.p2, "freq")
    print(d.temp, row.names = FALSE)
}
