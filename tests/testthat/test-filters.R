context("Filter functions")
test_that("test filter functions", {
  x <- filter_data(B2721)
  expect_is(x, "mappoly2.data")
  expect_is(x, "screened")
})
