context("training")

test_that("training successful", {
  data("parenthesis")
  preprocessed <- altum::preprocess(substr(parenthesis, 1, 100))
  expect_type(train_lstm(dat = preprocessed, vocabulary_size = 7), "list")
})