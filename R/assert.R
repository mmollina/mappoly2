is.mappoly2.data <- function (x)
  inherits(x, c("mappoly2.data"))

has.mappoly2.screened <- function(x){
  is.mappoly2.data(x) & inherits(x, c("screened"))
}

has.chromosome.info <- function(x){
  is.mappoly2.data(x) &
  !all(is.na(x$chrom)) & !is.null(x$chrom)
}

has.mappoly2.rf <- function(x){
  has.mappoly2.screened(x) & inherits(x, c("pairwise.rf"))
}

data.has.genome.info <- function(x){
  is.mappoly2.data(x) & !all(is.na(x$genome.pos))
}
is.mappoly2.sequence <- function (x)
  inherits(x, "mappoly2.sequence")

is.phased.sequence <- function(x){
  is.mappoly2.sequence(x) &&
    !is.null(x$phases[[1]]$p1) && !is.null(x$phases[[1]]$p2)
}


