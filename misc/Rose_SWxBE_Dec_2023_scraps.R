
layout(matrix(1:14, nrow = 2, byrow = T))
for(i in 1:7){
  id1 <- so[[i]]$maps$map.mds$order$locimap$locus[so[[i]]$maps$map.mds$order$locimap$confplotno]
  mappoly2:::plot_mappoly2_rf_matrix(so$data$pairwise.rf, ord = id1, fact = 5)
  id2 <- rownames(so[[i]]$maps$map.genome$order)
  mappoly2:::plot_mappoly2_rf_matrix(so$data$pairwise.rf, ord = id2, fact = 5)
}


layout(matrix(1:8, nrow = 2, byrow = T))
for(i in 1:7){
  id1 <- s[[i]]$maps$map.mds$order$locimap$locus[s[[i]]$maps$map.mds$order$locimap$confplotno]
  id2 <- rownames(s[[i]]$maps$map.genome$order)
  plot(match(id1,id2))
}
