context("predict")

test_that("prediction of next character", {
  example_model <- keras::load_model_hdf5("../../data-raw/example_model.hdf5") # model requires length of 30 char
  expect_error(predict_next_nucleotide())
  expect_error(predict_next_nucleotide(sequence = "", model = NULL))
  expect_s4_class(predict_next_nucleotide("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", example_model), "prediction")
})
