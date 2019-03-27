context("infer states")

test_that("infer states successful", {
  model_path <- "example_model.hdf5"
  preprocessed <- preprocess(strrep("A",80),
                             vocabulary = c("\n", "a", "c", "g", "t"))
  expect_type(getstates(model_path = model_path,
                        x = preprocessed$X,
                        maxlen = 30), "double")
})

