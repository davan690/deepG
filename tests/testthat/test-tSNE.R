context("tSNE")

test_that("Check Cell States", {
  expect_error(generateStatesFromFolder())
  expect_error(generateStatesFromFolder(""))
})

test_that("Check Cell States from FASTA files", {
  expect_error(extractCellFromStates())
  expect_error(extractCellFromStates(""))
})

test_that("Check Cell States from FASTA files", {
  expect_error(extractCellFromStates())
  expect_error(extractCellFromStates(""))
})
