library(dplyr)
require(visNetwork)
fl.names<-list.files("~/repos/official_repos/mappoly2/R", pattern = ".R", full.names = TRUE)
file.color <- RES <- NULL
for(j in 1:length(fl.names)){
  cat(paste(fl.names[j], "\n"))
  x<-readLines(con <- fl.names[j])
  doc.R <- !grepl(pattern = "#'", x = x)
  x <- x[doc.R]
  fc.start <- grep(pattern = "<- function", x = x)
  if(length(fc.start)==0) next()
  fc.start <- cbind(fc.start, c(fc.start[-1]-1, length(x)))
  A <- vector("list", nrow(fc.start))
  fc.name <- NULL
  for(i in 1:nrow(fc.start)){
    A[[i]]<-x[fc.start[i,1]:fc.start[i,2]]
    fc.name <- c(fc.name, strsplit(A[[i]][1], split = " ")[[1]][1])
  }
  names(A)<-fc.name
  RES <- c(RES, A)
  file.color <- c(file.color, rep(mappoly::mp_pallet3(length(fl.names))[j], length(A)))
}
names(file.color) <- all.fc.names <- names(RES)
df<-NULL
for(i in 1:length(RES)){
  for(j in 1:length(RES[[i]])){
    for(k in 1:length(all.fc.names)){
      if(names(RES[i]) != all.fc.names[k] & grepl(pattern = all.fc.names[k], RES[[i]][j]))
      {
        df <- rbind(df,data.frame(from = names(RES[i]), to = all.fc.names[k]))
      }
    }
  }
}
df<-unique(df)
df<-df[df$to != "##",]
head(df)
df
file.color[label]

nodes <- tibble(id = seq_along(label), label = label, color = file.color[label])
edges <- tibble(from = match(df[,1], label),
                to = match(df[,2], label))
visNetwork(nodes[nodes$label!="",], edges, width = "100%", height = "800px") %>%
  visEdges(arrows = "from")

