context("predict")

test_that("Prediction of next character", {
  #example.model <-
  #  keras::load_model_hdf5("example_model.hdf5") # model requires length of 30 char
  sequence <- strrep("A", 100)
  expect_error(predictNextNucleotide())
  expect_error(predictNextNucleotide(sequence = ""))
  expect_error(predictNextNucleotide(model = ""))
  #expect_type(predictNextNucleotide(sequence = sequence, model = example.model),"S4")
})

test_that("Prediction of replacement of n characters", {
  #example.model <-
  # keras::load_model_hdf5("example_model.hdf5")
  expect_error(replaceChar())
  expect_error(replaceChar(sequence = "", model = ""))
  #expect_type(replaceChar(sequence = sequence, model = example.model), "S4")
})