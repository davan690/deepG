context("preprocess")

test_that("correct vocabulary extraction", {
  expect_equal(get_vocabulary("ABC"), c("a", "b", "c"))
  expect_equal(get_vocabulary("CBA"), c("a", "b", "c"))
  expect_equal(get_vocabulary("AAA"), c("a"))
  expect_equal(get_vocabulary("012"), c("0", "1", "2"))
  expect_equal(get_vocabulary("()/"), c("(", ")", "/"))
  expect_equal(get_vocabulary(" a "), c(" ", "a"))
  expect_error(get_vocabulary(""))
  expect_error(get_vocabulary())
})

test_that("generating semi-redundant chunks", {
  expect_equivalent(lengths(preprocess("abcd", maxlen = 2))[1], 8)
  expect_equivalent(lengths(preprocess("abcd", maxlen = 2))[2], 4)
  expect_equivalent(preprocess("abcd", maxlen = 2)$Y, c(0,0,1,0))
})