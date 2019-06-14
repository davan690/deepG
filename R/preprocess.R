#' Returns the vocabulary from character string
#'
#' Use this function with a character string as input such as data(train)
#'
#' @param char a character string of text, length of one
#' @param verbose TRUE/FALSE
#' @export
get_vocabulary <- function(char, verbose = F) {
  library(dplyr)
  stopifnot(!is.null(char))
  stopifnot(nchar(char)>0)
  vocabulary <-  char %>%  stringr::str_to_lower() %>%
    stringr::str_c(collapse = "\n") %>%
    tokenizers::tokenize_characters(strip_non_alphanum = FALSE, simplify = TRUE) %>%
    unique() %>%
    sort()
  if (verbose) print(vocabulary)
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
#' @param maxlen length of semi-redundant sequences of maxlen characters
#' @param vocabulary char, should be sorted, if not set char vocabulary will be used
#' @param verbose TRUE/FALSE
#' @export
preprocess <- function(char, maxlen = 30, vocabulary, verbose=F) {
  require(dplyr)
  Check <- ArgumentCheck::newArgCheck()
  #* Add an error if maxlen <1
  if (maxlen < 1)
    ArgumentCheck::addError(
      msg = "'maxlen' must be >= 1",
      argcheck = Check
    )

  if (!missing(vocabulary)) {
    #* Add an error if vocabulary is <1
    if (length(vocabulary) < 1)
      ArgumentCheck::addError(
        msg = "'vocabulary' must be a character of length > 0",
        argcheck = Check
      )
  }
  #* Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(Check)

  # Load, collapse, and tokenize text
  text <- char %>%
    stringr::str_to_lower() %>%
    stringr::str_c(collapse = "\n") %>%
    tokenizers::tokenize_characters(strip_non_alphanum = FALSE, simplify = TRUE)
  if (verbose) print(sprintf("corpus length: %d", length(text)))

  if(missing(vocabulary)) {
    vocabulary <- text %>%
      unique() %>%
      sort()
  }
  if (verbose) print(sprintf("vocabulary size: %d", length(vocabulary)))

  # Cut the text in semi-redundant sequences of maxlen characters
  if (verbose) print("generation of semi-redundant sequences ...")
  dataset <- purrr::map(seq(1, length(text) - maxlen - 1, by = 1),
                        ~ list(sentece = text[.x:(.x + maxlen - 1)],
                               next_char = text[.x + maxlen]))
  dataset <- purrr::transpose(dataset)

  # Vectorization
  x <-
    array(0, dim = c(length(dataset$sentece), maxlen, length(vocabulary)))
  y <-
    array(0, dim = c(length(dataset$sentece), length(vocabulary)))
  if (verbose) print("vectorization ...")
  if (verbose) pb <-  txtProgressBar(min = 0,
                   max = length(dataset$sentece),
                   style = 3)
  for (i in 1:length(dataset$sentece)) {
  if (verbose) setTxtProgressBar(pb, i)
    x[i, ,] <- sapply(vocabulary, function(x) {
      as.integer(x == dataset$sentece[[i]])
    })
    y[i,] <- as.integer(vocabulary == dataset$next_char[[i]])
  }
  results <- list("X" = x, "Y" = y)
  return(results)
}

#' preprocesses and serializes genomes located in a folder, exports as hdf5
#'
#' @param directory input directory where .fasta files are located
#' @param out output directory wher hdf5 files will be written to
#' @param format either .txt or .fasta
#' @export
serialize_genomes <- function(directory, out, format = "fa"){
  require(xfun)
  require(rhdf5)
  library(Biostrings)
  directory <- normalize_path(directory)
  files <- list.files(path = directory, pattern = paste0("*.", format), full.names = TRUE)
  for (i in seq_along(files)) {
    file_name <- gsub(pattern = paste0("\\.", format, "$"), "", basename(files[i]))
   if (format == "fasta" | format == "fa") {
      fastaFile <- readDNAStringSet(files[i])
      seq <- paste(fastaFile) # not sure if this is working
      pre <-  preprocess(seq) # , ...)
   } else {
      stop("only .fa and .fasta files are supported")
    }
    #save
    filename <- paste0(out,"/",file_name, ".h5")
    h5createFile(filename) 
    h5createGroup(filename, "X")
    h5createGroup(filename, "Y")
    h5write(pre$X, filename, "X/data") 
    h5write(pre$Y, filename, "Y/data") 
    h5closeAll()
    }
}

#' wrapper of the preprocess() function called on the genomic contents of one
#' fasta file. Multiple entries are combined with newline characters.
#' @export
preprocess_fasta <- function(path, maxlen = 30,
                             vocabulary = c("\n", "a", "c", "g", "t"),
                             verbose = F){
  library(Biostrings)
  fastaFile <- readDNAStringSet(path)
  seq <- paste0(paste(fastaFile, collapse = "\n"), "\n")
  seq_processed <- preprocess(seq, maxlen = 30, vocabulary, verbose = F)
  return(seq_processed)
}

#' custom generator for fasta files
#'
#' @param directory input directory where .fasta files are located
#' @export
fasta_files_generator <- function(dir, format = "fasta") {
  library(xfun)
  fasta_files <- list.files(path = normalize_path(dir), 
                      pattern = paste0("*.", format), full.names = TRUE)
  next_file <- 0
  function() {
    # move to the next file (note the <<- assignment operator)
    next_file <<- next_file + 1
    # start from the beginning
    if (next_file > length(fasta_files))
      next_file <<- 1
    # determine the file name
    file <- fasta_files[[next_file]]
    preprocessed <- preprocess_fasta(file)
    # return the batch
    list(preprocessed$X, preprocessed$Y)
  }
}
