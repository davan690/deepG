#' Returns the vocabulary from character string
#'
#' Use this function with a character string
#'
#' @param char a character string of text with the length of one
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
    message("The vocabulary:", vocabulary)
  return(vocabulary)
}

#' Preprocess string to semi-redundant one-hot vector
#'
#' Outputs semi-redundant set of input character string.
#' Collapse, tokenize, and vectorize the character.
#' Use this function with a character string as input. For example, 
#' if the input text is ABCDEFGHI and the length(maxlen) is 5, the generating chunks would be:
#' X(1): ABCDE and Y(1): F;
#' X(2): BCDEF and Y(2): G;
#' X(3): CDEFG and Y(3): H;
#' X(4): DEFGH and Y(4): I
#' 
#' @param char a character input string of text with the length of one
#' @param labels a character string of same length as char with character as labels
#' @param maxlen length of the semi-redundant sequences
#' @param vocabulary char contains the vocabulary from the input char, it should be sorted
#' If no vocabulary exists, it is generated from the input char
#' @param verbose TRUE/FALSE
#' @example preprocessSemiRedudant("abcd",labels=NULL,maxlen=2,verbose=F)
#' @export
preprocessSemiRedundant <- function(char,
                                    labels = NULL,
                                    maxlen = 250,
                                    vocabulary,
                                    verbose = F) {
  require(dplyr)
  
  stopifnot(!is.null(char))
  stopifnot(nchar(char) > 0)
  stopifnot(maxlen > 0)
  
  # Load, collapse, and tokenize text ("ACGT" -> "a" "c" "g" "t")
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
    message("Finding the vocabulary ...")
  
  # Generating vocabulary from input char with the function getVocabulary()
  if (missing(vocabulary)) {
    vocabulary <- getVocabulary(char)
  }
  if(verbose)
    message("Vocabulary size:", length(vocabulary))
  # Cut the text in semi-redundant sequences of maxlen characters
  
  if (verbose)
    message("Generation of semi-redundant sequences ...")
  if (is.null(labels)) {
    dataset <- purrr::map(seq(1, length(text) - maxlen, by = 1),
                          ~ list(sentece = text[.x:(.x + maxlen - 1)],
                                 next_char = text[.x + maxlen]))
    dataset <- purrr::transpose(dataset)
    x <-
      array(0, dim = c(length(dataset$sentece), maxlen, length(vocabulary)))
    y <- array(0, dim = c(length(dataset$sentece), length(vocabulary)))
    # Vectorization
  if (verbose)
    message("Vectorization ...")
  if (verbose)
    pb <-  txtProgressBar(min = 0,
                          max = length(dataset$sentece),
                          style = 3)
    for (i in 1:length(dataset$sentece)) {
      if (verbose)
        setTxtProgressBar(pb, i)
      # generate one-hot encoding for one subset
      x[i, ,] <- sapply(vocabulary, function(x) {
        as.integer(x == dataset$sentece[[i]])
      })
      # target (next nucleotide in sequence)
      y[i,] <- as.integer(vocabulary == dataset$next_char[[i]])
    }
    } else {
      # use labels
      dataset <- purrr::map(seq(1, length(text) - maxlen, by = 1),
                            ~ list(sentece = text[.x:(.x + maxlen - 1)],
                                   label = text.labels[.x + maxlen]))
      dataset <- purrr::transpose(dataset)
      
      x <-
        array(0, dim = c(length(dataset$sentece), maxlen, length(vocabulary)))
      y <- array(0, dim = c(length(dataset$sentece), length(text.labels.vocabulary)))
      
      # Vectorization
      if (verbose)
        message("Vectorization ...")
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

#' Wrapper of the preprocessSemiRedundant()-function 
#' 
#' It called on the genomic contents of one
#' FASTA file. Multiple entries are combined with newline characters.
#' @param path path to the FASTA file
#' @param maxlen length of the semi-redundant sequences
#' @param vocabulary char contains the vocabulary
#' @param verbose TRUE/FALSE
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
    preprocessSemiRedundant(char = seq, maxlen = maxlen, vocabulary = vocabulary,
                            verbose = F) 
  } else {
    seq.processed <-
      preprocessSemiRedundant(char = seq, labels = seq.labels, maxlen = maxlen,
                              vocabulary = vocabulary, verbose = F) 
  }
  return(seq.processed)
}

#' Calculate the steps per epoch
#' 
#' Do one full preprocessing iteration to the FASTA file to figure out what the observed
#' steps_per_epoch value is.
#' @param dir 
#' @param batch.size
#' @param maxlen
#' @param format
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
      path = xfun::normalize_path(dir),
      pattern = paste0("*.", format),
      full.names = TRUE
    )
    for (file in fasta.files) {
      fasta.file <- Biostrings::readDNAStringSet(file)
      seq <- paste0(paste(fasta.file, collapse = "\n"), "\n")
      steps.per.epoch <-
        steps.per.epoch + ceiling((nchar(seq) - maxlen - 1) / batch.size)
    }
    return(steps.per.epoch)
  }

#' Custom generator for FASTA files
#' 
#' Produce chunks in size of batch.size
#' by iterating over the input files. If the input file is smaller than the
#' batch.size, batch size will be reduced to the maximal size. So the last batch
#' of a file is usually smaller.
#' See <https://github.com/bagasbgy/kerasgenerator/blob/master/R/timeseries_generator.R>
#' @param corpus.dir input directory where .fasta files are located
#' @param labels.dir
#' @param format
#' @param batch.size
#' @param maxlen
#' @param verbose TRUE/FALSE
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
    use_labels <<- TRUE
    preprocessed <- preprocessFasta(path = file, labels = label, maxlen = maxlen)
  } else {
    use_labels <<- FALSE
    preprocessed <- preprocessFasta(path = file, maxlen = maxlen)
  }
  if (verbose)
    message("Initializing ...")
  function() {
    # move to nexdt file if we cannot process another batch
    if (((batch.num) * batch.size) > nrow(preprocessed$X)) {
      if (verbose)
        message("Reached end of file.")
      # move to the next file
      next.file <<- next.file + 1
      # reset batch coordinates
      batch.num <<- 0
      batch.end <<- 0
      # at the end of the file list, start from the beginning
      if (next.file > length(fasta.files)) {
        if (verbose)
          message("Resetting to first file.")
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
        message("Reduce batch.size temporarily.")
      batch.end <- nrow(preprocessed$X)
      # reduced batch size
    } else {
      # regular batch.size
      batch.end <<- batch.start + batch.size  # 12
    }
    if (verbose)
      message(
        paste(
          "Generating batch number",
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
      preprocessed$X[batch.start:batch.end, , ]# dim should be (batch_size, maxlen, words)
    y.batch <-
      preprocessed$Y[batch.start:batch.end,]   # dim should be (batch_size, words)
    # return the file
    list(x.batch, y.batch)
  }
}
