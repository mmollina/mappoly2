#' Multipoint analysis using Hidden Markov Models
#'
#' @param void internal function
#' @keywords internal
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
mapping <- function(input.seq,
                    phase.conf = "all",
                    rf = NULL,
                    error = 0.0,
                    verbose = FALSE,
                    tol = 10e-4,
                    ret_H0 = FALSE)
{
  assert_that(is.mappoly2.sequence(input.seq))
  mrk.id <- rownames(input.seq$phases[[1]]$p1)
  g <- input.seq$data$geno.dose[mrk.id, ]
  g[is.na(g)] <- -1
  if(is.null(rf))
    rf <- rep(0.01, nrow(g) - 1)
  assert_that(length(rf) == nrow(g) - 1)
  if(all(phase.conf == "all"))
    phase.conf <- 1:length(input.seq$phases)

  assert_that(all(phase.conf%in%1:length(input.seq$phases)),
              msg = "invalid phases specified in 'phase.conf'")

  output.seq <- input.seq

  cat("Multi-locus map estimation\n")
  cat("   Number of phase configurations: ", length(phase.conf), "\n")
  if (detect_info_par(input.seq) == "both" || ret_H0){
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
                                 err = error,
                                 verbose = verbose,
                                 detailed_verbose = FALSE,
                                 tol = tol,
                                 ret_H0 = ret_H0)
      output.seq$phases[[i]]$loglike <- w[[1]]
      output.seq$phases[[i]]$rf <- w[[2]]
      output.seq$phases[[i]]$error <- error
    }
    cat("Done with map estimation\n")

    return(sort_phase(output.seq))
  } else if(detect_info_par(input.seq) == "p1"){
    id <- which(input.seq$data$ploidy.p2 == input.seq$data$dosage.p2[mrk.id])
    g[id, ] <- g[id, ] - input.seq$data$ploidy.p2/2
    for(i in phase.conf){
      cat("   Conf.", i,":")
      w <- est_hmm_map_biallelic_single(PH = input.seq$phases[[i]]$p1,
                                        G = g,
                                        rf = rf,
                                        err = error,
                                        verbose = verbose,
                                        detailed_verbose = FALSE,
                                        tol = tol,
                                        ret_H0 = ret_H0)
      output.seq$phases[[i]]$loglike <- w[[1]]
      output.seq$phases[[i]]$rf <- w[[2]]
      output.seq$phases[[i]]$error <- error
    }
    cat("Done with map estimation\n")
    return(sort_phase(output.seq))
  } else if(detect_info_par(input.seq) == "p2") {
    id <- which(input.seq$data$ploidy.p1 == input.seq$data$dosage.p1[mrk.id])
    g[id, ] <- g[id, ] - input.seq$data$ploidy.p1/2
    for(i in phase.conf){
      cat("   Conf.", i,":")
      w <- est_hmm_map_biallelic_single(PH = input.seq$phases[[i]]$p2,
                                        G = g,
                                        rf = rf,
                                        err = error,
                                        verbose = verbose,
                                        detailed_verbose = FALSE,
                                        tol = tol,
                                        ret_H0 = ret_H0)
      output.seq$phases[[i]]$loglike <- w[[1]]
      output.seq$phases[[i]]$rf <- w[[2]]
      output.seq$phases[[i]]$error <- error
    }
    cat("Done with map estimation\n")
    return(sort_phase(output.seq))
  } else {
    stop("it should not get here")
  }
}
#' Efficiently phases unprocessed markers in a given phased
#' genetic map, circumventing the need for complete
#' HMM-based map recomputation.
#'
#' @param void internal function
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
augment_phased_map <- function(input.seq,
                               input.twopt,
                               thresh.LOD.ph = 5,
                               thresh.LOD.rf = 5,
                               thresh.rf = 0.5,
                               max.phases = 5,
                               thresh.LOD.ph.to.insert = 10,
                               thresh.rf.to.insert = NULL,
                               tol = 10e-4,
                               verbose = TRUE){
  assert_that(is.haplotype.sequence(input.seq))
  if(all(input.seq$mrk.names%in%rownames(input.seq$phases[[1]]$p1))){
    message("All markers are phased. Retutning original sequence.")
    return(input.seq)
  }
  assert_that(is.mappoly2.twopt(input.twopt))
  M <- rf_list_to_matrix(input.twopt,
                         thresh.LOD.ph = thresh.LOD.ph,
                         thresh.LOD.rf = thresh.LOD.rf,
                         thresh.rf = thresh.rf,
                         shared.alleles = TRUE)
  assert_that(matrix_contain_data_seq(M,input.seq))
  mrk.pos <- rownames(input.seq$phases[[1]]$p1) # positioned markers
  mrk.id <- setdiff(input.seq$mrk.names, mrk.pos) # markers to be positioned
  ## two-point phasing parent 1
  dose.vec <- input.seq$data$dosage.p1[mrk.id]
  InitPh1 <- input.seq$phases[[1]]$p1
  S1 <- M$Sh.p1[mrk.id, mrk.pos]
  L1 <- mappoly2:::phasing_one(mrk.id, dose.vec, S1, InitPh1, verbose)
  ## two-point phasing parent 2
  dose.vec <- input.seq$data$dosage.p2[mrk.id]
  InitPh2 <- input.seq$phases[[1]]$p2
  S2 <- M$Sh.p2[mrk.id, mrk.pos]
  L2 <- mappoly2:::phasing_one(mrk.id, dose.vec, S2, InitPh2, verbose)
  ## Selecting phase configurations
  n.conf <- sapply(L1, nrow) * sapply(L2, nrow)
  if(verbose){
    cat("Distribution of phase configurations.\n")
    temp <- as.data.frame(table(n.conf))
    colnames(temp) <- c("n. phase. conf.", "frequency")
    print(temp)
  }
  mrk.sel <- which(n.conf <= max.phases)
  if(length(mrk.sel) == 0)
    stop("No markers were selected for 'max.phases' = ", max.phases,
         "\n'max.phases' should be at least ", min(n.conf))
  L1 <- L1[n.conf <= max.phases]
  L2 <- L2[n.conf <= max.phases]
  mrk.id <- mrk.id[n.conf <= max.phases]
  pedigree <- matrix(rep(c(1,
                           2,
                           input.seq$data$ploidy.p1,
                           input.seq$data$ploidy.p2, 1),
                         input.seq$data$n.ind),
                     nrow = input.seq$data$n.ind,
                     byrow = TRUE)
  flanking <- mappoly2:::find_flanking_markers(input.seq$mrk.names, mrk.pos, mrk.id)
  phasing_results <- vector("list", length(flanking))
  names(phasing_results) <- names(flanking)
  if(verbose) pb <- txtProgressBar(min = 0, max = length(L1), style = 3)
  for(i in 1:length(L1)){
    G <- input.seq$data$geno.dose[mrk.id[i], ,drop = TRUE]
    G[is.na(G)] <- -1
    u <- match(unlist(flanking[[mrk.id[i]]]), mrk.pos)
    if(is.na(u)[1]){ # Marker inserted at the beginning of the linkage group
      homolog_prob <- as.matrix(input.seq$phases[[1]]$haploprob[,c(na.omit(u), na.omit(u)+1)+2])
      idx <- c(1,0,2)
    } else if(is.na(u)[2]){ # Marker inserted at the end of the linkage group
      homolog_prob <- as.matrix(input.seq$phases[[1]]$haploprob[,c(na.omit(u)-1, na.omit(u))+2])
      idx <- c(0,2,1)
    } else { # Marker inserted in the middle of the linkage group
      homolog_prob <- as.matrix(input.seq$phases[[1]]$haploprob[,u+2])
      idx <- c(0,1,2)
    }
    w2<-w1<-NULL
    z<-vector("list", nrow(L1[[i]]) * nrow(L2[[i]]))
    count <- 1
    for(j in 1:nrow(L1[[i]])){
      for(k in 1:nrow(L2[[i]])){
        PH <- list(L1[[i]][j,], L2[[i]][k,])
        z[[count]]<-mappoly2:::est_hmm_map_biallelic_insert_marker(PH,
                                                                   G,
                                                                   pedigree,
                                                                   homolog_prob,
                                                                   rf = c(0.01,0.01),
                                                                   idx,
                                                                   verbose = FALSE,
                                                                   detailed_verbose = FALSE,
                                                                   tol = tol,
                                                                   ret_H0 = FALSE)
        w1 <- rbind(w1, L1[[i]][j,])
        w2 <- rbind(w2, L2[[i]][k,])
        count <- count + 1
      }
    }
    x <- sapply(z, function(x) x[[1]])
    x <- max(x) - x
    id <- order(x)
    phasing_results[[mrk.id[i]]] <- list(loglike = x[id],
                                         rf.vec = t(sapply(z[id],
                                                           function(x) x[[2]])),
                                         phases = list(p1 = w1[id,,drop=FALSE],
                                                       p2 = w2[id,,drop=FALSE]))
    if(verbose) setTxtProgressBar(pb, i)

  }
  if(verbose) close(pb)
  selected.list<-phasing_results[sapply(phasing_results,
                                        function(x) length(x$loglike)==1 ||
                                          x$loglike[2] > thresh.LOD.ph.to.insert)]
  if(is.null(thresh.rf.to.insert))
    thresh.rf.to.insert <- max(input.seq$phases[[1]]$rf)
  if(thresh.rf.to.insert < 0 || thresh.rf.to.insert >= 0.5)
    stop("'thresh.rf.to.insert' parameter must be between 0 and 0.5")
  selected.list <- phasing_results[sapply(phasing_results, function(x) max(x$rf.vec[1,]) <= thresh.rf.to.insert)]
  for(i in names(selected.list)){
    pos <- mappoly2:::find_flanking_markers(input.seq$mrk.names,
                                            rownames(input.seq$phases[[1]]$p1),
                                            i)
    if(length(unlist(pos)) == 0) next()
    cur.mrk <- rownames(input.seq$phases[[1]]$p1)
    if(is.na(pos[[1]]$preceding))# beginning
    {
      input.seq$phases[[1]]$p1 <- rbind(selected.list[[i]]$phases$p1[1,], input.seq$phases[[1]]$p1)
      input.seq$phases[[1]]$p2 <- rbind(selected.list[[i]]$phases$p2[1,], input.seq$phases[[1]]$p2)
      rownames(input.seq$phases[[1]]$p1) <- rownames(input.seq$phases[[1]]$p2) <- c(i, cur.mrk)
    }
    else if (is.na(pos[[1]]$succeeding)){ #end
      input.seq$phases[[1]]$p1 <- rbind(input.seq$phases[[1]]$p1, selected.list[[i]]$phases$p1[1,])
      input.seq$phases[[1]]$p2 <- rbind(input.seq$phases[[1]]$p2, selected.list[[i]]$phases$p2[1,])
      rownames(input.seq$phases[[1]]$p1) <- rownames(input.seq$phases[[1]]$p2) <- c(cur.mrk, i)
    }
    else {
      preceding <- cur.mrk[1:match(pos[[1]]$preceding, cur.mrk)]
      succeeding <- cur.mrk[(match(pos[[1]]$succeeding, cur.mrk)): length(cur.mrk)]
      input.seq$phases[[1]]$p1 <- rbind(input.seq$phases[[1]]$p1[preceding,],
                                        selected.list[[i]]$phases$p1[1,],
                                        input.seq$phases[[1]]$p1[succeeding,])
      input.seq$phases[[1]]$p2 <- rbind(input.seq$phases[[1]]$p2[preceding,],
                                        selected.list[[i]]$phases$p2[1,],
                                        input.seq$phases[[1]]$p2[succeeding,])
      rownames(input.seq$phases[[1]]$p1) <- rownames(input.seq$phases[[1]]$p2) <- c(preceding, i, succeeding)
    }
  }
  return(input.seq)
}

#'This function merges two genetic maps built from markers that are exclusively
#'informative in isolated parents, facilitating unified analysis and visualization
#'of distinct genetic data.
#'
#' @param void internal function
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
merge_single_parent_maps <- function(input.seq.all,
                                     input.seq.p1,
                                     input.seq.p2,
                                     input.twopt){
  temp.mrk.id <- character(length(input.seq.all$mrk.names))
  temp.mrk.id[match(rownames(input.seq.p1$phases[[1]]$p1), input.seq.all$mrk.names)] <- rownames(input.seq.p1$phases[[1]]$p1)
  temp.mrk.id[match(rownames(input.seq.p2$phases[[1]]$p2), input.seq.all$mrk.names)] <- rownames(input.seq.p2$phases[[1]]$p2)
  temp.mrk.id <- temp.mrk.id[temp.mrk.id != ""]
  ph.p2 <- ph.p1 <- matrix(0, length(temp.mrk.id), input.seq.p1$data$ploidy.p1, dimnames = list(temp.mrk.id, NULL))
  ph.p1[rownames(input.seq.p1$phases[[1]]$p1), ] <- input.seq.p1$phases[[1]]$p1
  ph.p2[rownames(input.seq.p1$phases[[1]]$p2), ] <- input.seq.p1$phases[[1]]$p2
  ph.p1[rownames(input.seq.p2$phases[[1]]$p1), ] <- input.seq.p2$phases[[1]]$p1
  ph.p2[rownames(input.seq.p2$phases[[1]]$p2), ] <- input.seq.p2$phases[[1]]$p2
  input.seq.all$phases <- list(list(p1 = ph.p1,
                                    p2 = ph.p2,
                                    loglike = NULL,
                                    rf = NULL,
                                    error = NULL,
                                    haploprob = NULL))
  return(input.seq.all)
}

#' prepare maps for plot
#' @param void internal function to be documented
#' @keywords internal
prepare_map <- function(input.seq){
  assert_that(mappoly2:::is.mapped.sequence(input.seq))
  ## Gathering marker positions
  map <- cumsum(imf_h(c(0, input.seq$phases[[1]]$rf)))
  names(map) <- rownames(input.seq$phases[[1]]$p1)

  ## Gathering phases
  ph.p1 <- input.seq$phases[[1]]$p1
  ph.p2 <- input.seq$phases[[1]]$p2
  colnames(ph.p1) <- paste0("p1.", 1:input.seq$data$ploidy.p1)
  colnames(ph.p2) <- paste0("p2.", 1:input.seq$data$ploidy.p2)
  if(is.null(input.seq$data$ref))
  {
    ph.p1[ph.p1 == 1] <- ph.p2[ph.p2 == 1] <- "A"
    ph.p1[ph.p1 == 0] <- ph.p2[ph.p2 == 0] <- "B"
  } else {
    for(i in names(map)){
      ph.p1[i, ph.p1[i,] == 1] <- input.seq$data$alt[i]
      ph.p1[i, ph.p1[i,] == 0] <- input.seq$data$ref[i]
      ph.p2[i, ph.p2[i,] == 1] <- input.seq$data$alt[i]
      ph.p2[i, ph.p2[i,] == 0] <- input.seq$data$ref[i]
    }
  }
  d.p1 <- input.seq$data$dosage.p1[names(map)]
  d.p2 <- input.seq$data$dosage.p2[names(map)]
  list(ploidy.p1 = input.seq$data$ploidy.p1,
       ploidy.p2 = input.seq$data$ploidy.p2,
       name.p1 = input.seq$data$name.p1,
       name.p2 = input.seq$data$name.p2,
       map = map,
       ph.p1 = ph.p1,
       ph.p2 = ph.p2,
       d.p1 = d.p1,
       d.p2 = d.p2)
}

#' @importFrom grDevices rgb
#' @importFrom graphics rect
#' @export
plot_map <- function(x, left.lim = 0, right.lim = Inf,
                     phase = TRUE, mrk.names = FALSE,
                     plot.dose = TRUE, homolog.names.adj = 3,
                     cex = 1, xlim = NULL, main = "",...) {
  old.p1ar <- par(no.readonly = TRUE)
  on.exit(par(old.p1ar))
  map.info <- prepare_map(x)
  if(any(map.info$ph.p1 == "B")){
    var.col <- c(A = "black", B = "darkgray")
  } else {
    var.col <- c(A = "#008000", T = "#FF0000", C = "#0000FF", G = "#FFFF00")
  }
  ploidy <- max(c(map.info$ploidy.p1, map.info$ploidy.p2))
  x <- map.info$map
  lab <- names(x)
  zy <- seq(0, 0.6, by = 0.12)
  zy.p1 <- zy[1:map.info$ploidy.p1] +1.8 + (0.3 * ((map.info$ploidy.p2/2)-1))
  zy.p2 <- zy[1:map.info$ploidy.p2] + 1.1
  pp <- map.info$ph.p1
  pq <- map.info$ph.p2
  d.p1 <- map.info$d.p1
  d.p2 <- map.info$d.p2
  x1 <- abs(left.lim - x)
  x2 <- abs(right.lim - x)
  id.left <- which(x1 == min(x1))[1]
  id.right <- rev(which(x2 == min(x2)))[1]
  par(mai = c(1,0.15,0,0), mar = c(4.5,homolog.names.adj,1,2))
  curx <- x[id.left:id.right]
  layout(mat  = matrix(c(2,4,1,3), ncol = 2), heights = c(10, 1), widths = c(1, 10))
  #layout(mat  = matrix(c(4,2,3, 1), ncol = 2), heights = c(2, 10), widths = c(1, 10))
  if(is.null(xlim)){
    xlim <- range(curx)
  }
  max.y <- 4.0
  plot(x = curx,
       y = rep(.5,length(curx)),
       type = "n" ,
       ylim = c(.25, max.y),
       axes = FALSE,
       xlab = "Distance (cM)",
       ylab = "",
       xlim = xlim)
  lines(c(x[id.left], x[id.right]), c(.5, .5), lwd = 15, col = "gray")
  points(x = curx,
         y = rep(.5,length(curx)),
         xlab = "", ylab = "",
         pch = "|", cex = 1.5,
         ylim = c(0,2))
  axis(side = 1)
  ####Parent 2#####
  x1 <- seq(x[id.left], x[id.right], length.out = length(curx))
  x.control <- diff(x1[1:2])/2
  if(length(x1) < 150)
    x.control <- x.control * .8
  if(length(x1) < 100)
    x.control <- x.control * .8
  if(length(x1) < 75)
    x.control <- x.control * .8
  if(length(x1) < 50)
    x.control <- x.control * .8
  if(length(x1) < 25)
    x.control <- x.control * .8
  for(i in 1:map.info$ploidy.p2)
  {
    lines(range(x1), c(zy.p2[i], zy.p2[i]), lwd = 12, col = "gray")
    y1 <- rep(zy.p2[i], length(curx))
    pal <- var.col[pq[id.left:id.right,i]]
    rect(xleft = x1 - x.control,
         ybottom = y1 -.035,
         xright = x1 + x.control,
         ytop = y1 +.035,
         col = pal,
         border = NA)
  }
  #connecting allelic variants to markers
  for(i in 1:length(x1))
    lines(c(curx[i], x1[i]), c(0.575, zy.p2[1]-.05), lwd = 0.2)
  ####
  if(plot.dose){
    y <- zy.p2[map.info$ploidy.p2]+0.1 - ((map.info$ploidy.p2/2 - 1)*0.005)+d.p2[id.left:id.right]/20
    y.l <- zy.p2[map.info$ploidy.p2]+0.1 - ((map.info$ploidy.p2/2 - 1)*0.005)+ c(0:map.info$ploidy.p2)[id.left:id.right]/20
    #text(x = min(x1) - 1, y = mean(y.l), "doses", srt = 90)
    for(i in 1:length(y.l)){
      text(x = min(x1)-1,y=y.l[i],i-1, cex = .7)
      lines(range(x1), c(y.l[i], y.l[i]), lwd = .5, col = "gray")
    }

    points(x = x1,
           y = y,
           col = "darkgray",
           #col = d.col[as.character(d.p2[id.left:id.right])],
           pch = 19, cex = .7)
  }
  ####Parent 1#####
  for(i in 1:map.info$ploidy.p1)
  {
    lines(range(x1), c(zy.p1[i], zy.p1[i]), lwd = 12, col = "gray")
    y1 <- rep(zy.p1[i], length(curx))
    pal <- var.col[pp[id.left:id.right,i]]
    rect(xleft = x1 - x.control,
         ybottom = y1 -.035,
         xright = x1 + x.control,
         ytop = y1 +.035,
         col = pal,
         border = NA)
  }
  ####
  if(plot.dose){
    y <- zy.p1[map.info$ploidy.p1]+0.1 -((map.info$ploidy.p1/2 - 1)*0.005)+d.p1[id.left:id.right]/20
    y.l <- zy.p1[map.info$ploidy.p1]+0.1 -((map.info$ploidy.p1/2 - 1)*0.005)+c(0:map.info$ploidy.p1)[id.left:id.right]/20
    #text(x = min(x1) - 1, y = mean(y.l), "doses", srt = 90)
    for(i in 1:length(y.l)){
      text(x = min(x1)-1,y=y.l[i],i-1, cex = .7)
      lines(range(x1), c(y.l[i], y.l[i]), lwd = .5, col = "gray")
    }
    points(x = x1,
           y = y,
           col = "darkgray",
           #col = d.col[as.character(d.p1[id.left:id.right])],
           pch = 19, cex = .7)
  }
  if(mrk.names)
    text(x = x1,
         y = rep(max(y)+.1, length(x1)),
         labels = names(curx),
         srt = 90, adj = 0, cex = cex *.6)
  par(mar = c(4.5,1,1,0), xpd = TRUE)
  plot(x = 0,
       y = 0,
       type = "n" ,
       axes = FALSE,
       ylab = "",
       xlab = "",
       ylim = c(.25, max.y))

  mtext(text = main, side = 2, at = mean(c(zy.p2, zy.p2)), line = -1, font = 4, cex = cex , adj = c(0,0))
  mtext(text = map.info$name.p2, side = 4, at = mean(zy.p2), line = -1, font = 4)
  for(i in 1:map.info$ploidy.p2)
    mtext(colnames(map.info$ph.p2)[i], line = 1, at = zy.p2[i], side = 4, las = 2, cex = 0.7 * cex)
  mtext(text = map.info$name.p1, side = 4, at = mean(zy.p1), line = -1, font = 4)
  for(i in 1:map.info$ploidy.p1)
    mtext(colnames(map.info$ph.p1)[i],  line = 1, at = zy.p1[i], side = 4, las = 2, cex = 0.7 * cex)
  par(mar = c(0,0,0,0), xpd = FALSE)
  plot(x = curx,
       y = rep(.5,length(curx)),
       type = "n" ,
       axes = FALSE,
       xlab = "",
       ylab = "")
  if(any(map.info$ph.p1 == "B")){
    legend("topleft", legend = c("A", "B"),
           fill  = c(var.col), #title = "Variants",
           box.lty = 0, bg = "transparent", ncol = 6)
  } else {
    legend("topleft", legend = c("A", "T", "C", "G", "-"),
           fill  = c(var.col, "white"),# title = "Nucleotides",
           box.lty = 0, bg = "transparent", ncol = 6)
  }
}

#' Physical versus genetic distance
#'
#' This function plots scatterplot(s) of physical distance (in Mbp) versus the genetic
#' distance (in cM). Map(s) should be passed as a single object or a list of objects
#' of class \code{mappoly.map}.
#'
#' @param  seq.list A list or a single object of class \code{mappoly.map}
#'
#' @param phase.config A vector containing which phase configuration should be
#'  plotted. If \code{'best'} (default), plots the configuration
#'  with the highest likelihood for all elements in \code{'seq.list'}
#'
#' @param same.ch.lg Logical. If \code{TRUE} displays only the scatterplots between the
#'   chromosomes and linkage groups with the same number. Default is \code{FALSE}.
#'
#' @param alpha transparency factor for SNPs points
#'
#' @param size size of the SNP points
#'
#' @examples
#'   plot_genome_vs_map(solcap.mds.map, same.ch.lg = TRUE)
#'   plot_genome_vs_map(solcap.mds.map, same.ch.lg = FALSE,
#'                      alpha = 1, size = 1/2)
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_.
#'     \doi{10.1534/g3.119.400378}
#'
#' @export plot_genome_vs_map
plot_genome_vs_map <- function(seq.list,
                               same.ch.lg = FALSE,
                               alpha = 1/5,
                               size = 3){
  if(mappoly2:::is.mappoly2.sequence(seq.list)){
    assert_that(mappoly2:::is.mapped.sequence(seq.list))
    seq.list <- list(seq.list)
  }
  assert_that(is.list(seq.list))
  if (any(!sapply(seq.list, mappoly2:::is.mapped.sequence)))
    stop("All elemnts in 'seq.list' should be mapped sequences.")
  if(any(sapply(seq.list, function(x) is.null(x$data$genome.pos))))
    stop("All elemnts in 'seq.list' should have the genomic position of the markers")
  geno.vs.map <- NULL
  for(i in 1:length(seq.list)){
    LG <- genomic.pos <- map.pos <- NULL
    mrk.names <- rownames(seq.list[[i]]$phases[[1]]$p1)
    geno.vs.map <- rbind(geno.vs.map,
                         data.frame(mrk.names = mrk.names,
                                    map.pos = cumsum(imf_h(c(0, seq.list[[i]]$phases[[1]]$rf))),
                                    genomic.pos = seq.list[[i]]$data$genome.pos[mrk.names]/1e6,
                                    LG = as.factor(i),
                                    chr = as.factor(seq.list[[i]]$data$chrom[mrk.names])))
  }
  geno.vs.map$chr <- factor(geno.vs.map$chr, levels = sort(levels(geno.vs.map$chr)))
  if(same.ch.lg){
    p <- ggplot2::ggplot(geno.vs.map, ggplot2::aes(genomic.pos, map.pos)) +
      ggplot2::geom_point(alpha = alpha, ggplot2::aes(colour = LG), size = size) +
      ggplot2::facet_wrap(~LG, nrow = floor(sqrt(length(seq.list)))) +
      ggplot2::labs(subtitle = "Linkage group", x = "Genome position (Mbp)", y = "Map position (cM)") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     legend.position = "none", plot.subtitle = ggplot2::element_text(hjust = 0.5))

  } else {
    p <- ggplot2::ggplot(geno.vs.map, ggplot2::aes(genomic.pos, map.pos)) +
      ggplot2::geom_point(alpha = alpha, ggplot2::aes(colour = LG), size = size) +
      ggplot2::facet_grid(LG~chr) +
      ggplot2::labs(x = "Genome position (Mbp)", y = "Map position (cM)") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "none")
  }
  p
}
