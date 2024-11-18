# Sample test for `read_geno_csv`
test_that("read_geno_csv processes real CSV file correctly", {
  # Locate the CSV file from the package's extdata directory
  tempfl <- list.files(system.file("extdata", package = "mappoly2"), full.names = TRUE)

  # Run the function
  result <- read_geno_csv(
    file.in = tempfl,
    ploidy.p1 = 4,
    name.p1 = "I195",
    name.p2 = "F1.85.209"
  )

  # Ensure the result is not empty and is of the correct class
  expect_s3_class(result, "mappoly2.data")
  expect_true(result$n.mrk > 0) # Check non-zero markers
  expect_true(result$n.ind > 0) # Check non-zero individuals

  # Ensure parental information is correct
  expect_equal(result$name.p1, "I195")
  expect_equal(result$name.p2, "F1.85.209")
})
test_that("read_vcf processes real VCF file correctly", {
  # Prepare the VCF file
  fl <- "https://github.com/mmollina/MAPpoly_vignettes/raw/master/data/sweet_sample_ch3.vcf.gz"
  tempfl <- tempfile(pattern = 'chr3_', fileext = '.vcf.gz')
  download.file(fl, destfile = tempfl)

  # Run the function
  result <- read_vcf(
    file.in = tempfl,
    ploidy.p1 = 6,
    name.p1 = "PARENT1",
    name.p2 = "PARENT2",
    min.av.depth = 20 # Minimum average depth to retain markers
  )

  # Ensure the result is not empty and is of the correct class
  expect_s3_class(result, "mappoly2.data")
  expect_true(result$n.mrk > 0) # Check non-zero markers
  expect_true(result$n.ind > 0) # Check non-zero individuals

  # Clean up the temporary file
  unlink(tempfl)
})

