context("infer-states")

test_that("Check Cell States", {
  expect_error(getStates())
  expect_error(getStates(model.path = ""))
  expect_error(getStates(x = ""))
  expect_error(getStates(model.path = "", x = ""))
})

test_that("Check Cell States from FASTA files", {
  expect_error(getStatesFromFasta())
})
