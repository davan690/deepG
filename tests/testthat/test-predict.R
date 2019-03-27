context("predict")

test_that("prediction of next character", {
  example_model <- keras::load_model_hdf5("example_model.hdf5") # model requires length of 30 char
  expect_error(predict_next_nucleotide())
  expect_error(predict_next_nucleotide(sequence = "", model = NULL))
  expect_s4_class(predict_next_nucleotide(strrep("A",30),
                                          example_model), "prediction")
})

test_that("prediction of replacement of n characters", {
  example_model <- keras::load_model_hdf5("example_model.hdf5") # model requires length of 30 char
  expect_error(replace_char())
  expect_error(replace_char(sequence = "", model = NULL))
  expect_type(replace_char(strrep("A",30),
                                          example_model), "character")
})
