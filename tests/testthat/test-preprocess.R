context("preprocess")

test_that("correct vocabulary extraction", {
  expect_equal(getVocabulary("ABC"), c("a", "b", "c"))
  expect_equal(getVocabulary("CBA"), c("a", "b", "c"))
  expect_equal(getVocabulary("AAA"), c("a"))
  expect_equal(getVocabulary("012"), c("0", "1", "2"))
  expect_equal(getVocabulary("()/"), c("(", ")", "/"))
  expect_equal(getVocabulary(" a "), c(" ", "a"))
  expect_error(getVocabulary(""))
  expect_error(getVocabulary())
})

test_that("generating semi-redundant chunks", {
  expect_equivalent(lengths(preprocessSemiRedundant(char = "abcd", maxlen = 2))[1], 8)
  expect_equivalent(lengths(preprocessSemiRedundant(char = "abcd", maxlen = 2))[2], 4)
  expect_equivalent(preprocessSemiRedundant(char = "abcd", maxlen = 2)$Y, c(0,0,1,0))
})
