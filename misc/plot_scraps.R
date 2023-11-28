plot.mappoly2 <- function(x,
                          what = c("raw",
                                   "screened",
                                   "initiated",
                                   "pairwise",
                                   "pairwise_mds",
                                   "pairwise_genome",
                                   "group"),
                          type = c("rf", "lod"),
                          ord = NULL,
                          rem = NULL,
                          main.text = NULL,
                          index = FALSE,
                          fact = 1, ...)
{
  oldpar <- par()
  on.exit(par(oldpar))
  what <- match.arg(what)
  if (what == "raw"){
    assert_that(inherits(x, "data"))
    plot_data(x, text = "Raw data", col = "darkred")
  }
  else if(what == "screened"){
    assert_that(inherits(x, "screened"))
    plot_data(x,
              text = "All data - screened",
              col =  "#07607e",
              mrk.id = x$screened.data$mrk.names,
              ind.id = x$screened.data$ind.names)
  }
  else if(what == "initiated"){
    assert_that(inherits(x, "initiated"))
    plot_data(x,
              text = "Initial sequence - screened",
              col =  "#07607e",
              mrk.id = x$initial.sequence,
              ind.id = x$screened.data$ind.names)
  }
  else if(what == "pairwise"){
    assert_that(inherits(x, "pairwise"))
    plot_rf_matrix(x$pairwise,
                   type = type,
                   ord = ord,
                   rem = rem,
                   main.text = main.text,
                   index = index,
                   fact = fact)
  } else if(what == "group"){
    plot_group(x$linkage.groups)
  } else if(what == "pairwise_mds"){
    assert_that(inherits(x, "pairwise"))
    idx <- sapply(x$working.sequences, function(x) !is.null(x$order$mds))
    assert_that(any(idx))
    if(sum(idx) == 2)
      op <- par(mfrow = c(1, 2), pty = "s")
    else
      op <- par(mfrow = c(ceiling(sqrt(sum(idx))),
                          ceiling(sqrt(sum(idx)))),
                pty = "s")
    on.exit(par(op))
    for(i in which(idx))
      mappoly2:::plot_rf_matrix(x$pairwise,
                                type = type,
                                ord = x$working.sequences[[i]]$order$mds$info$locimap$locus,
                                main.text = paste0(names(x$working.sequences)[i], "-MDS"),
                                fact = fact)
  } else if(what == "pairwise_genome"){
    assert_that(inherits(x, "pairwise"))
    idx <- sapply(x$working.sequences, function(x) !is.null(x$order$genome))
    assert_that(any(idx))
    if(sum(idx) == 2)
      op <- par(mfrow = c(1, 2), pty = "s")
    else
      op <- par(mfrow = c(ceiling(sqrt(sum(idx))),
                          ceiling(sqrt(sum(idx)))),
                pty = "s")
    for(i in which(idx))
      mappoly2:::plot_rf_matrix(x$pairwise,
                                type = type,
                                ord = rownames(x$working.sequences[[i]]$order$genome$info),
                                main.text = paste0(names(x$working.sequences)[i], "-genome"),
                                fact = fact)
  }
}
