is.mappoly2.data <- function (x)
  inherits(x, c("mappoly2.data"))
is.mappol2.screened <- function(x)
  "screened.data"%in%names(x)
is.mappoly2.map <- function (x)
  inherits(x, "mappoly2.map")
is.mappoly2.sequence <- function (x)
  inherits(x, "mappoly2.sequence")
is.mappoly2.twopt <- function (x)
  inherits(x, "mappoly2.twopt")
is.mappoly2.pcmap <- function (x)
  inherits(x, "mappoly2.pcmap")
is.mappoly2.pcmap3d <- function (x)
  inherits(x, "mappoly2.pcmap3d")
is.mappoly2.group <- function (x)
  inherits(x, "mappoly2.group")
is.mappoly2.chitest.seq <- function (x)
  inherits(x, "mappoly2.chitest.seq")
is.mappoly2.geno.ord <- function (x)
  inherits(x, "mappoly2.geno.ord")
is.mappoly2.rf.matrix <- function (x)
  inherits(x, "mappoly2.rf.matrix")
is.pairwise.sequence <- function(x){
  is.mappoly2.sequence(x) &&
    !is.null(x$pairwise) &&
    all(colnames(x$pairwise$rec.mat) == x$mrk.names)
}
is.phased.sequence <- function(x){
  is.mappoly2.sequence(x) &&
  !is.null(x$phases[[1]]$p1) && !is.null(x$phases[[1]]$p2)
}
has.chromosome.info <- function(x){
  !all(is.na(x$data$chrom)) & !is.null(x$data$chrom)
}
is.grouped.sequence <- function(x){
  is.mappoly2.sequence(x) &&
    !is.null(x$linkage.groups)
}
is.mapped.sequence <- function(x){
  is.mappoly2.sequence(x) &&
  !is.null(x$phases[[1]]$rf)
}
is.haplotype.sequence <- function(x){
  is.mappoly2.sequence(x) &&
    !is.null(x$phases[[1]]$haploprob)
}
matrix_contain_data_seq <- function(M,x){
  all(x$data$mrk.names%in%rownames(M$Sh.p1))
}

