library(dplyr)
require(visNetwork)
fl.names<-list.files("~/repos/official_repos/mappoly2/R", pattern = ".R", full.names = TRUE)
RES <- NULL
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
}
all.fc.names <- names(RES)
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
df<-df[!apply(df, 1, function(x) any(grepl("is.mappoly2", x))),]
df<-df[!apply(df, 1, function(x) any(grepl("print.mappoly2", x))),]
df<-df[!apply(df, 1, function(x) any(grepl("plot.mappoly2", x))),]
df<-df[df$to != "##",]
# create a dataset:

#pal<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
pal<-c("#4363d8", "#f58231", '#ffe119','#e6194b','#3cb44b')
label <- unique(as.character(as.matrix(df)))
palette<-rep(pal[1], length(label))
matches <- grep(paste("read", "table_to",sep="|"),
                label)
palette[matches] <- pal[2]
matches <- grep(paste("filter",sep="|"),
                label)
palette[matches] <- pal[3]
matches <- grep(paste("rf", "pairwise",sep="|"),
                label)
palette[matches] <- pal[4]
matches <- grep(paste("group",sep="|"),
                label)
palette[matches] <- pal[5]

df <- df[label!="",]
palette <- palette[label!=""]
label <- label[label!=""]


nodes <- tibble(id = seq_along(label), label = label, color = palette)
edges <- tibble(from = match(df[,1], label),
                to = match(df[,2], label))
visNetwork(nodes, edges, width = "100%", height = "800px") %>%
  visEdges(arrows = "from")

