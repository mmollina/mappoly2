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

is.phased.sequence <- function(x, lg, type, parent, phase.type){
  assert_that(length(lg) == 1)
  assert_that(length(type) == 1)
  assert_that(length(parent) == 1)
  assert_that(length(phase.type) == 1)
  !is.null(x$maps[[lg]][[type]][[parent]][[phase.type]][[1]]$p1)
}

is.mds.ordered <- function(x, lg){
  assert_that(length(lg) == 1)
  !is.null(x$maps[[lg]][["mds"]]$order)
}

is.mapped.sequence <- function(x, lg, type, parent){
  assert_that(length(lg) == 1)
  assert_that(length(type) == 1)
  assert_that(length(parent) == 1)
  !is.null(x$maps[[lg]][[type]][[parent]]$hmm.phase[[1]]$loglike)
}

is.haplotype.sequence <- function(x, lg, type, parent){
  assert_that(length(lg) == 1)
  assert_that(length(type) == 1)
  assert_that(length(parent) == 1)
  !is.null(x$maps[[lg]][[type]][[parent]]$hmm.phase[[1]]$haploprob)
}

has.adequate.ploidy <- function(x){
  assertthat::assert_that((x %% 2) == 0, msg = "At least one of the parents has an odd ploidy level.")
}
