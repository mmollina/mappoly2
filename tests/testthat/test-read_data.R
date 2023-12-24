context("Read data")
test_that("read data from CSV file correctly", {
  fpath <- system.file("extdata", "B2721.csv", package="mappoly2")
  dat <- read_geno_csv(file.in  = fpath,
                       ploidy.p1 = 4,
                       ploidy.p2 = 4,
                       name.p1 = "Atlantic",
                       name.p2 = "B1829-5")
  expect_identical(dat, B2721)
  testthat_print(dat)
  x <- sapply(dat, length)
  y <- unique(x[c("mrk.names",
                  "dosage.p1",
                  "dosage.p2",
                  "chrom",
                  "genome.pos",
                  "ref",
                  "alt")])
  expect_equal(y, dat$n.mrk)
  expect_true(x["ind.names"] == dat$n.ind)
  expect_equal(length(y), 1)
  expect_true(x["geno.dose"]/y == x["ind.names"])
})

