#' Returns the vocabulary from character string
#'
#' Use this function with a character string as input such as data(train)
#'
#' @param char a character string of text, length of one
#' @param verbose TRUE/FALSE
#' @export
getVocabulary <- function(char, verbose = F) {
  library(dplyr)
  stopifnot(!is.null(char))
  stopifnot(nchar(char) > 0)
  vocabulary <-  char %>%  stringr::str_to_lower() %>%
    stringr::str_c(collapse = "\n") %>%
    tokenizers::tokenize_characters(strip_non_alphanum = FALSE, simplify = TRUE) %>%
    unique() %>%
    sort()
  if (verbose)
    print(vocabulary)
  return(vocabulary)
}

#' Preprocess string to semi redundant one-hot vector
#'
#' Outputs semi-redundant set of input string
#' Preprocess a character input
#' Collapse and tokenize and vectorize character
#' Use this function with a character string as input such as data(train)
#' if the input text ist ABCDEFGHIJKLM and the maxlen is set to 5, the chunks would be
#' X(1): ABCDE Y(1):F
#' X(2): BCDEF Y(2):G
#' X(3): CDEFG Y(3):H
#' X(4): DEFGH (4):I
#' ...
#'
#' @param char a character string of text, length of one
#' @param labels a character string of same length as char with character as labels
#' @param maxlen length of semi-redundant sequences of maxlen characters
#' @param vocabulary char, should be sorted, if not set char vocabulary will be used
#' @param verbose TRUE/FALSE
#' @export
preprocessSemiRedundant <- function(char,
                                    labels = NULL,
                                    maxlen = 250,
                                    vocabulary,
                                    verbose = F) {
  require(dplyr)
  Check <- ArgumentCheck::newArgCheck()
  #* Add an error if maxlen <1
  if (maxlen < 1)
    ArgumentCheck::addError(msg = "'maxlen' must be >= 1",
                            argcheck = Check)
  
  if (!missing(vocabulary)) {
    #* Add an error if vocabulary is <1
    if (length(vocabulary) < 1)
      ArgumentCheck::addError(msg = "'vocabulary' must be a character of length > 0",
                              argcheck = Check)
  }
  #* Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(Check)
  
  # Load, collapse, and tokenize text
  text <- char %>%
    stringr::str_to_lower() %>%
    stringr::str_c(collapse = "\n") %>%
    tokenizers::tokenize_characters(strip_non_alphanum = FALSE, simplify = TRUE)
  
  if (!is.null(labels)) {
    text.labels <- labels %>%
      stringr::str_c(collapse = "\n") %>%
      tokenizers::tokenize_characters(strip_non_alphanum = FALSE, simplify = TRUE)
    text.labels.vocabulary <- text.labels %>%
      unique() %>%
      sort()
  }
  
  if (verbose)
    print(sprintf("corpus length: %d", xt))
  
  if (missing(vocabulary)) {
    vocabulary <- text %>%
      unique() %>%
      sort()
  }
  if (verbose)
    print(sprintf("vocabulary size: %d", length(vocabulary)))
  
  # Cut the text in semi-redundant sequences of maxlen characters
  if (verbose)
    print("generation of semi-redundant sequences ...")
  if (is.null(labels)) {
    dataset <- purrr::map(seq(1, length(text) - maxlen - 1, by = 1),
                          ~ list(sentece = text[.x:(.x + maxlen - 1)],
                                 next_char = text[.x + maxlen]))
    dataset <- purrr::transpose(dataset)
    x <-
      array(0, dim = c(length(dataset$sentece), maxlen, length(vocabulary)))
    y <- array(0, dim = c(length(dataset$sentece), length(vocabulary)))
    # Vectorization
    if (verbose)
      print("vectorization ...")
    if (verbose)
      pb <-  txtProgressBar(min = 0,
                            max = length(dataset$sentece),
                            style = 3)
    for (i in 1:length(dataset$sentece)) {
      if (verbose)
        setTxtProgressBar(pb, i)
      x[i, ,] <- sapply(vocabulary, function(x) {
        as.integer(x == dataset$sentece[[i]])
      })
      y[i,] <- as.integer(vocabulary == dataset$next_char[[i]])
    }
    } else {
      # use labels
      dataset <- purrr::map(seq(1, length(text) - maxlen - 1, by = 1),
                            ~ list(sentece = text[.x:(.x + maxlen - 1)],
                                   label = text.labels[.x + maxlen]))
      dataset <- purrr::transpose(dataset)
      
      x <-
        array(0, dim = c(length(dataset$sentece), maxlen, length(vocabulary)))
      y <- array(0, dim = c(length(dataset$sentece), length(text.labels.vocabulary)))
      
      # Vectorization
      if (verbose)
        print("vectorization ...")
      if (verbose)
        pb <-  txtProgressBar(min = 0,
                              max = length(dataset$sentece),
                              style = 3)
      for (i in 1:length(dataset$sentece)) {
        if (verbose)
          setTxtProgressBar(pb, i)
        x[i, ,] <- sapply(vocabulary, function(x) {
          as.integer(x == dataset$sentece[[i]])
        })
        y[i,] <- as.integer(text.labels.vocabulary == dataset$label[[i]])
      }
    }

  results <- list("X" = x, "Y" = y)
  return(results)
}

#' wrapper of the preprocessSemiRedundant() function called on the genomic contents of one
#' fasta file. Multiple entries are combined with newline characters.
#' @export
preprocessFasta <- function(path,
                            labels = NULL,
                            maxlen = 250,
                            vocabulary = c("\n", "a", "c", "g", "t"),
                            verbose = F) {
  # process labels
  if (!is.null(labels)){
    fasta.file.labels <- Biostrings::readDNAStringSet(labels)
    seq.labels <- paste0(paste(fasta.file.labels, collapse = "\n"), "\n")
  } 
  
  # process corpus
  fasta.file <- Biostrings::readDNAStringSet(path)
  seq <- paste0(paste(fasta.file, collapse = "\n"), "\n")
  
  if (is.null(labels)){
  seq.processed <-
    preprocessSemiRedundant(char = seq, maxlen = maxlen, vocabulary, verbose = F) 
  } else {
    seq.processed <-
      preprocessSemiRedundant(char = seq, labels = seq.labels, maxlen = maxlen, vocabulary, verbose = F) 
  }
  return(seq.processed)
}

#' do one full preprocessing iteration to figure out what the observed
#' steps_per_epoch value is
#' @export
calculateStepsPerEpoch <-
  function(dir,
           batch.size = 256,
           maxlen = 250,
           format = "fasta") {
    library(xfun)
    library(Biostrings)
    steps.per.epoch <- 0
    fasta.files <- list.files(
      path = normalize_path(dir),
      pattern = paste0("*.", format),
      full.names = TRUE
    )
    for (file in fasta.files) {
      fasta.file <- readDNAStringSet(file)
      seq <- paste0(paste(fasta.file, collapse = "\n"), "\n")
      steps.per.epoch <-
        steps.per.epoch + ceiling((nchar(seq) - maxlen - 1) / batch.size)
    }
    return(steps.per.epoch)
  }

#' custom generator for fasta files, will produce chunks in size of batch.size
#' by iterating over the input files. If the input file is smaller than the
#' batch.size, batch size will be reduced to the maximal size. So the last batch
#' of a file is usually smaller.
#'
#'see https://github.com/bagasbgy/kerasgenerator/blob/master/R/timeseries_generator.R
#' @param directory input directory where .fasta files are located
#' @export
fastaFileGenerator <- function(corpus.dir,
                               labels.dir,
                               format = "fasta",
                               batch.size = 512,
                               maxlen = 250,
                               verbose = F) {
  fasta.files <- list.files(
    path = xfun::normalize_path(corpus.dir),
    pattern = paste0("*.", format),
    full.names = TRUE
  )
  
  if (!missing(labels.dir)){
    label.files <- paste0(labels.dir, gsub(pattern = paste0("\\.",format,"$"), "", basename(fasta.files)),".txt")
  }
  
  next.file <- 1
  batch.num <- 0
  batch.end <- 0
  # pre-load the first file
  file <- fasta.files[[1]]
  if (!missing(labels.dir)) {
    label <- label.files[[1]]
    use_labels <- TRUE
    preprocessed <- preprocessFasta(path = file, labels = label, maxlen = maxlen)
  } else {
    preprocessed <- preprocessFasta(path = file, maxlen = maxlen)
  }
  if (verbose)
    message("initializing")
  function() {
    # move to nexdt file if we cannot process another batch
    if (((batch.num) * batch.size) > nrow(preprocessed$X)) {
      if (verbose)
        message("reached end of file")
      # move to the next file
      next.file <<- next.file + 1
      # reset batch coordinates
      batch.num <<- 0
      batch.end <<- 0
      # at the end of the file list, start from the beginning
      if (next.file > length(fasta.files)) {
        if (verbose)
          message("resetting to first file")
        # reset file number
        next.file <<- 1
        # reset batch coordinates
        batch.num <<- 0
        batch.end <<- 0
      }
      # read in the new file
      file <<- fasta.files[[next.file]]
      if (use_labels) {
        preprocessed <- preprocessFasta(path = file, labels = label, maxlen = maxlen)
      } else {
        preprocessed <- preprocessFasta(path = file, maxlen = maxlen)
      }
    }
    # proceccing a batch
    batch.num <<- batch.num + 1
    batch.start <<- batch.end + 1
    # check if full batch size can be processed, if not temporaly reduce
    # batch_size to maximum possible
    if ((batch.start + batch.size) > nrow(preprocessed$X)) {
      if (verbose)
        message("reduce batch.size temporarily")
      batch.end <- nrow(preprocessed$X)
      # reduced batch size
    } else {
      # regular batch.size
      batch.end <<- batch.start + batch.size  # 12
    }
    if (verbose)
      message(
        paste(
          "generating bach number",
          batch.num,
          batch.start,
          "-",
          batch.end,
          "of file",
          file,
          "index",
          next.file
        )
      )
    x.batch <-
      preprocessed$X[batch.start:batch.end, , ]# dim shoiuld be (batch_size, length, words)
    y.batch <-
      preprocessed$Y[batch.start:batch.end,] #  # dim should be (batch_size, words)
    # return the file
    list(x.batch, y.batch)
  }
}
