genes <- c("MAPK1", "MAPK3")

test_that("genes not mapped returns null", {
  expect_null(compartmentData(c("madeupgene1", "madeupgene2")))
  expect_null(compartmentData(genes, organism = "Yeast"))
})

test_that("returns data frame",{
  expect_equal(class(compartmentData(genes)), "data.frame")
  expect_equal(class(compartmentData(genes, bkgd = genes)), "data.frame")
  expect_equal(class(compartmentData(genes, trafficking = T)), "data.frame")
  expect_equal(class(compartmentData(genes, subAnnots = T)), "data.frame")
})

test_that("returns error if empty vector",{
  expect_error(compartmentData(""))
})
