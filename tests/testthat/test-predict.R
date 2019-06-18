context("predict")


#test_that("prediction of next character", {
#  example.model <-
#    keras::load_model_hdf5("example_model.hdf5") # model requires length of 30 char
#  expect_error(predictNextNucleotide())
#  expect_error(predictNextNucleotide(sequence = "", model = NULL))
#  expect_s4_class(predictNextNucleotide(strrep("A", 30),
#                                        example.model), "prediction")
#})

test_that("prediction of replacement of n characters", {
  example.model <-
    keras::load_model_hdf5("example_model.hdf5") # model requires length of 30 char
  expect_error(replaceChar())
  expect_error(replaceChar(sequence = "", model = NULL))
  expect_type(replaceChar(strrep("A", 30),
                          example.model), "character")
})
