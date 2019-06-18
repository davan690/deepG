context("training")

test_that("training successful", {
  data("parenthesis")
  maxlen <- 30
  preprocessed <- altum::preprocessSemiRedundant(substr(parenthesis, 1, 100), 
                                                 maxlen = maxlen)
  expect_type(trainNetwork(dataset = preprocessed, 
                           vocabulary.size = 7, 
                           batch.size = 10,
                           maxlen = maxlen,
                           layers.lstm = 2,
                           layer.size = 10,), "list")
})
