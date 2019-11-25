context("predict")

test_that("Prediction of next character", {
  
  skip_if_no_keras()
  
  example.model <-
    keras::load_model_hdf5("example_model.hdf5") # model requires length of 30 char
  sequence <- strrep("A", 100)
  
  expect_error(predictNextNucleotide())
  expect_error(predictNextNucleotide(sequence = ""))
  expect_error(predictNextNucleotide(model = ""))
  
  predicted_NextNucleotide <- predictNextNucleotide(sequence = sequence, model = example.model, verbose = T)
  expect_message(predictNextNucleotide(sequence = sequence, model = example.model, verbose = T))
  expect_silent(predictNextNucleotide(sequence = sequence, model = example.model))
  expect_s4_class(predicted_NextNucleotide, "prediction")
  expect_type(predicted_NextNucleotide@next_char, "character")
  expect_type(predicted_NextNucleotide@probability, "double")
  expect_type(predicted_NextNucleotide@index, "integer")
  expect_type(predicted_NextNucleotide@alternative_probability, "double")
  expect_type(predicted_NextNucleotide@solution, "character")
  
})

test_that("Prediction of replacement of n characters", {
  
  skip_if_no_keras()
  
  example.model <-
    keras::load_model_hdf5("example_model.hdf5")
  sequence <- strrep("A", 100)
  
  expect_error(replaceChar())
  expect_error(replaceChar(sequence = "", model = ""))
  expect_type(replaceChar(sequence = sequence, model = example.model), "character")
  expect_equivalent(nchar(replaceChar(sequence = sequence, model = example.model)), 100)
  
  expect_silent(replaceChar(sequence = sequence, model = example.model))
})
