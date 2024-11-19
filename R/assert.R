# Function to check if an object is of class "mappoly2.data"
is.mappoly2.data <- function(x) {
  inherits(x, "mappoly2.data")
}

# Function to check if an object has been screened
has.mappoly2.screened <- function(x) {
  is.mappoly2.data(x) && length(x$screened.data) > 0
}

# Function to check if chromosome information is present
has.chromosome.info <- function(x) {
  is.mappoly2.data(x) &&
    !all(is.na(x$chrom)) && !is.null(x$chrom)
}

# Function to check if recombination frequency data is present
has.mappoly2.rf <- function(x) {
  has.mappoly2.screened(x) && inherits(x, "pairwise.rf")
}

# Function to check if genome position information is present
data.has.genome.info <- function(x) {
  is.mappoly2.data(x) && !all(is.na(x$genome.pos))
}

# Function to check if an object is a "mappoly2.sequence"
is.mappoly2.sequence <- function(x) {
  inherits(x, "mappoly2.sequence")
}

# Function to check if a sequence is phased
is.phased.sequence <- function(x) {
  is.mappoly2.sequence(x) &&
    !is.null(x$phases[[1]]$p1) &&
    !is.null(x$phases[[1]]$p2)
}

# Function to check if a sequence is phased based on specific parameters
is.phased.sequence <- function(x, lg, type, parent, phase.type) {
  assertthat::assert_that(length(lg) == 1)
  assertthat::assert_that(length(type) == 1)
  assertthat::assert_that(length(parent) == 1)
  assertthat::assert_that(length(phase.type) == 1)
  !is.null(x$maps[[lg]][[type]][[parent]][[phase.type]][[1]]$p1)
}

# Function to check if a sequence is MDS ordered
is.mds.ordered <- function(x, lg) {
  assertthat::assert_that(length(lg) == 1)
  !is.null(x$maps[[lg]][["mds"]]$order)
}

# Function to check if a sequence is mapped
is.mapped.sequence <- function(x, lg, type, parent) {
  assertthat::assert_that(length(lg) == 1)
  assertthat::assert_that(length(type) == 1)
  assertthat::assert_that(length(parent) == 1)
  !is.null(x$maps[[lg]][[type]][[parent]]$hmm.phase[[1]]$loglike)
}

# Function to check if a sequence is a haplotype sequence
is.haplotype.sequence <- function(x, lg, type, parent) {
  assertthat::assert_that(length(lg) == 1)
  assertthat::assert_that(length(type) == 1)
  assertthat::assert_that(length(parent) == 1)
  !is.null(x$maps[[lg]][[type]][[parent]]$hmm.phase[[1]]$haploprob)
}

# Function to check if ploidy levels are adequate (even)
has.adequate.ploidy <- function(x) {
  assertthat::assert_that(
    (x %% 2) == 0,
    msg = "At least one of the parents has an odd ploidy level."
  )
}
