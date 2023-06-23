twopt_phasing <- function(mrk.id, # vector of marker names of size m
                          dose, # vector of integers of size m
                          ploidy, # even number representing ploidy level (2,4,6)
                          S, # matrix of integers varying from 0 to ploidy, with dimensions (m x m)
                          Ph){
  H <- list(Ph[1,,drop = FALSE])
  idx <- 1
  for(i in 2:length(mrk.id)){
    x <- dose[i]
    if(x == ploidy){
      for(j in 1:length(H)){
        H[[j]] <- rbind(H[[j]], rep(1, ploidy))
      }
      idx <- c(idx, i)
    }
    else if(x == 0){
      for(j in 1:length(H)){
        H[[j]] <- rbind(H[[j]], rep(0, ploidy))
      }
      idx <- c(idx, i)
    }
    else {
      d <-  S[i, 1:(i-1)]
      vtemp <- vector("list", length(H))
      Hres <- NULL
      for(j in 1:length(H)){
        vtemp[[j]] <- mappoly2:::find_valid_permutations(H = H[[j]],
                                                         d = d[idx],
                                                         x = x)
        if(nrow(vtemp[[j]]) == 0) next()
        Htemp <- vector("list", nrow(vtemp[[j]]))
        for(k in 1:nrow(vtemp[[j]])){
          Htemp[[k]] <- rbind(H[[j]], vtemp[[j]][k,,drop = FALSE])
        }
        Hres <- c(Hres, mappoly2:::filter_matrices(Htemp))
      }
      if(length(Hres) > 1){
        H <- mappoly2:::filter_matrices(Hres)
        idx <- c(idx, i)
      } else if(length(Hres) == 1) {
        H <- Hres
        idx <- c(idx, i)
      }
    }
  }
  for(i in 1:length(H))
    rownames(H[[i]]) <- mrk.id[idx]
  return(H)
}

compare_phases <- function(X,Y){
  pl <- ncol(X)
  if(!identical(dim(X),dim(Y))) stop()
  I <- gtools::permutations(pl,pl,1:pl)
  r <- numeric(nrow(I))
  for(i in 1:nrow(I)){
    r[i] <- cor(as.numeric(Y[,I[i,]]), as.numeric(X))
  }
  return(max(abs(r)))
}

