context("preprocess")

test_that("NAs are not dropped from the data", {
  expect_equal(get_vocabulary("ABC"), c("a", "b", "c"))
  expect_equal(get_vocabulary("CBA"), c("a", "b", "c"))
  expect_equal(get_vocabulary("AAA"), c("a"))
  expect_equal(get_vocabulary("012"), c("0", "1", "2"))
  expect_equal(get_vocabulary("()/"), c("(", ")", "/"))
  expect_equal(get_vocabulary(" a "), c(" ", "a"))
  expect_error(get_vocabulary(""))
  expect_error(get_vocabulary())
})
