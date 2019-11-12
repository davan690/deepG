context("preprocess")

test_that("Correct vocabulary extraction", {
  
  expect_equal(getVocabulary("ABC"), c("a", "b", "c"))
  expect_equal(getVocabulary("CBA"), c("a", "b", "c"))
  expect_equal(getVocabulary("AAA"), c("a"))
  expect_equal(getVocabulary("012"), c("0", "1", "2"))
  expect_equal(getVocabulary("()/"), c("(", ")", "/"))
  expect_equal(getVocabulary(" a "), c(" ", "a"))
  expect_equal(getVocabulary("abc"), c("a","b","c"))
  expect_equal(getVocabulary(123),c("1","2","3"))
  expect_equal(getVocabulary(c("A","B","C",1)),c("\n","1","a","b","c"))
  expect_equal(getVocabulary("\n"),c("\n"))
  
  expect_error(getVocabulary(""))
  expect_error(getVocabulary())
  
  expect_is(getVocabulary("abc"),"character")
  
  expect_message(getVocabulary("abc", verbose = T))
  expect_silent(getVocabulary("abc"))
})

test_that("Generating semi-redundant chunks", {
  
  expect_is(preprocessSemiRedundant(char = "abcd", maxlen = 2),"list")
  expect_is(preprocessSemiRedundant(char = "abcd", maxlen = 2)$X,"array")
  expect_is(preprocessSemiRedundant(char = "abcd", maxlen = 2)$Y,"matrix")
  
  expect_equivalent(lengths(preprocessSemiRedundant(char = "abcd", maxlen = 2))[1], 16)
  expect_equivalent(lengths(preprocessSemiRedundant(char = "abcd", maxlen = 2, vocabulary = c("a","b","c","d")))[1], 16)
  expect_equivalent(lengths(preprocessSemiRedundant(char = "abcd", maxlen = 2))[2], 8)
  expect_equivalent(preprocessSemiRedundant(char = "abcd", maxlen = 2)$Y, matrix(c(0,0,1,0,0,0,0,1), byrow = TRUE, nrow = 2))           
  expect_equivalent(length(preprocessSemiRedundant(char="abcd", maxlen = 2)),2)
  
  expect_error(preprocessSemiRedundant(char = "abcd", maxlen = ""))
  expect_error(preprocessSemiRedundant(char = "abcd", maxlen = 0))
  expect_error(preprocessSemiRedundant(char = "abcd", vocabulary = ""))
  expect_error(preprocessSemiRedundant(char = "abcd", vocabulary = 0))
  
  expect_message(preprocessSemiRedundant(char = "abcd", maxlen = 2, verbose = T))
  expect_silent(preprocessSemiRedundant(char = "abcd", maxlen = 2))

  expect_type(preprocessSemiRedundant(char = "abcd", maxlen = 2)$X, "double")
  expect_type(preprocessSemiRedundant(char = "abcd", maxlen = 2)$Y, "double")
})

test_that("Generating semi-redundant chunks from Fasta files", {
  
  file <- file.path("fasta/a.fasta")
  
  expect_is(preprocessFasta(file),"list")
  expect_is(preprocessFasta(file)$X,"array")
  expect_is(preprocessFasta(file)$Y,"matrix")
  
  expect_equivalent(lengths(preprocessFasta(file))[2], 39672)
  expect_equivalent(length(preprocessFasta(file)),2)
  
  expect_error(preprocessFasta())
  expect_error(preprocessFasta(""))
  
  expect_silent(preprocessFasta(file))
  
  expect_type(preprocessFasta(file)$X, "double")
  expect_type(preprocessFasta(file)$Y, "double")
})

test_that("Checking the generator for the Fasta files", {
  
  testpath <- file.path("fasta/")
  batch.size = 80
  maxlen = 50
  words = 6

  gen <- fastaFileGenerator(corpus.dir =  testpath, batch.size = batch.size, maxlen = maxlen, seqStart = "", showWarnings = FALSE)
  
  expect_equivalent(dim(gen()[[1]])[1], batch.size)
  expect_equivalent(dim(gen()[[1]])[2], maxlen)
  expect_equivalent(dim(gen()[[1]])[3], words)
  expect_equivalent(dim(gen()[[2]])[1], batch.size)
  expect_equivalent(dim(gen()[[2]])[2], words)
  expect_equivalent(length(gen()),2)
  
  expect_error(fastaFileGenerator())
  expect_error(fastaFileGenerator(""))
  expect_error(fastaFileGenerator(testpath, batch.size = batch.size, maxlen = maxlen, seqStart = "l",
                                  vocabulary = c("x","y","z")))
  
  expect_is(fastaFileGenerator(testpath, batch.size = batch.size, maxlen = maxlen, seqStart = "", showWarnings = FALSE), "function")
  expect_is(gen(), "list")
  expect_is(gen()[[1]], "array")
  expect_is(gen()[[2]], "matrix")
  
  expect_message(fastaFileGenerator(testpath, batch.size = batch.size, maxlen = maxlen, 
                                    seqStart = "", showWarnings = FALSE, verbose = T))
  expect_silent(fastaFileGenerator(testpath, batch.size = batch.size, maxlen = maxlen, seqStart = "", showWarnings = FALSE))
  
  expect_warning(fastaFileGenerator(testpath, batch.size = batch.size, maxlen = maxlen, seqStart = "", seqEnd= "",
                                              withinFile = "", vocabulary = c("x","y","z", "a", "b", "c"), showWarnings = TRUE))
  
  
  expect_type(gen()[[1]], "integer")
  expect_type(gen()[[2]], "integer")
})
