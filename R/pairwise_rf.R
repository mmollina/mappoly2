
#' @rdname rf_list_to_matrix
#' @export
plot_mappoly2_rf_matrix <- function(x, type = c("rf", "lod"), ord = NULL, rem = NULL,
                                    main.text = NULL, index = FALSE, fact = 1, ...){
  type <- match.arg(type)
  if(is.mappoly2.sequence(ord))
    ord <- ord$mrk.names
  if(type  ==  "rf"){
    w <- x$rec.mat
    if(!is.null(ord))
    {
      w <- w[ord,ord]
    }
    if(!(is.null(rem) || sum(colnames(x$rec.mat)%in%rem)  ==  0))
    {
      o <- which(colnames(x$rec.mat)%in%rem)
      w <- w[-o,-o]
    }
    if(fact > 1)
      w <- aggregate_matrix(w, fact)
    if(is.null(main.text))
      main.text <- "Recombination fraction matrix"
    col.range  <-
      na.omit(rev(fields::tim.colors())[1:(ceiling(128 * max(x$rec.mat, na.rm = TRUE)) + 1)])
    brks <- NULL
  } else if(type  ==  "lod")
  {
    w <- x$lod.mat
    if(!is.null(ord))
    {
      w <- w[ord,ord]
    }
    if(!(is.null(rem) || sum(colnames(x$rec.mat)%in%rem)  ==  0))
    {
      o <- which(colnames(x$rec.mat)%in%rem)
      w <- w[-o,-o]
    }
    if(fact > 1)
      w <- aggregate_matrix(w, fact)
    w[w < 1e-4] <- 1e-4
    w <- log10(w)
    if(is.null(main.text))
      main.text <- "log(LOD) Score matrix"
    col.range <- na.omit(fields::tim.colors()[1:(ceiling(128 * max(x$lod.mat, na.rm = TRUE)) + 1)])
    col.range <- col.range[ceiling(seq(1, length(col.range), length.out = 10))]
    brks <- seq(min(w, na.rm = TRUE), max(w, na.rm = TRUE), length.out = 11)
    brks <- round(exp(brks/log10(exp(1))),1)
  } else stop("Invalid matrix type.")

  fields::image.plot(
    w,
    col = col.range,
    lab.breaks = brks,
    main = main.text,
    useRaster = FALSE,
    axes = FALSE
  )
  if(ncol(w) < 100)
    ft <- .7
  else
    ft <- 100/ncol(w)
  if(index)
    text(x = seq(0,1, length.out = ncol(w)), y = seq(0,1, length.out = ncol(w)),
         labels = colnames(w), cex = ft)
}


#' Select rf and lod based on thresholds
#'
#' @param void internal function to be documented
#' @keywords internal
select_rf <- function(x, thresh.LOD.ph, thresh.LOD.rf, thresh.rf, shared.alleles = FALSE)
{
  if(any(is.na(x)))
  {
    if(shared.alleles){return(c(NA,NA,NA,NA,NA))} else return(c(NA,NA,NA))
  }
  if((nrow(x)  ==  1 || abs(x[2 , 1]) >= thresh.LOD.ph) &&
     abs(x[1, 3]) >= thresh.LOD.rf &&
     abs(x[1, 2]) <= thresh.rf)
  {
    if(shared.alleles){
      y <- strsplit(rownames(x), "-")
      return(c(x[1,2:3], x[2,1], as.numeric(y[[1]][1]), as.numeric(y[[1]][2])))
    } else {
      return(c(x[1,2:3], x[2,1]))
    }
  }
  else{
    {
      if(shared.alleles){return(c(NA,NA,NA,NA,NA))} else return(c(NA,NA,NA))
    }
  }
}


#' Aggregate matrix cells (lower the resolution by a factor)
#'
#' @param void internal function to be documented
#' @keywords internal
#' @export
aggregate_matrix <- function(M, fact){
  id <- seq(1,ncol(M), by = fact)
  id <- cbind(id, c(id[-1]-1, ncol(M)))
  R <- matrix(NA, nrow(id), nrow(id))
  for(i in 1:(nrow(id)-1)){
    for(j in (i+1):nrow(id)){
      R[j,i] <-  R[i,j] <- mean(M[id[i,1]:id[i,2], id[j,1]:id[j,2]], na.rm = TRUE)
    }
  }
  R
}
