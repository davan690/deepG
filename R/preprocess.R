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

#' do one full preprocessing iteration to figure out what the observed
#' steps_per_epoch value is
calculate_steps_per_epoch <- function(dir, batch_size = 12){
  library(xfun)
  steps_per_epoch <- 0
  fasta_files <- list.files(path = normalize_path(dir), 
                            pattern = paste0("*.", format), full.names = TRUE)
  for (file in fasta_files){
    preprocessed <- preprocess_fasta(file)
    steps_per_epoch <- steps_per_epoch +  ceiling(nrow(preprocessed$X)/batch_size)
  }
  return(steps_per_epoch)
}

#' custom generator for fasta files, will produce chunks in size of batch_size
#' by iterating over the input files. If the input file is smaller than the
#' batch_size, batch size will be reduced to the maximal size. So the last batch
#' of a file is usually smaller. 
#'
#'see https://github.com/bagasbgy/kerasgenerator/blob/master/R/timeseries_generator.R
#' @param directory input directory where .fasta files are located
#' @export
fasta_files_generator <- function(dir, format = "fasta", 
                                  batch_size = 12,
                                  verbose = F) {
  library(xfun)
  fasta_files <- list.files(path = normalize_path(dir), 
                      pattern = paste0("*.", format), full.names = TRUE)
  next_file <- 1
  batch_num <- 0
  batch_end <- 0
  # pre-load the first file
  file <- fasta_files[[1]]
  preprocessed <- preprocess_fasta(file)
  if (verbose) message("initializing")
  function() {
    # move to nexdt file if we cannot process another batch
    if (((batch_num) * batch_size) > nrow(preprocessed$X)){
      if (verbose) message("reached end of file")
      # move to the next file
      next_file <<- next_file + 1 
      # reset batch coordinates
      batch_num <<- 0
      batch_end <<- 0
      # at the end of the file list, start from the beginning
      if (next_file > length(fasta_files)) {
        if (verbose) message("resetting to first file")
        # reset file number
        next_file <<- 1
        # reset batch coordinates
        batch_num <<- 0
        batch_end <<- 0
      }
      # read in the new file
      file <<- fasta_files[[next_file]]
      preprocessed <- preprocess_fasta(file)
    }
    # proceccing a batch
    batch_num <<- batch_num + 1
    batch_start <<- batch_end + 1
    # check if full batch size can be processed, if not temporaly reduce
    # batch_size to maximum possible
    if ((batch_start + batch_size) > nrow(preprocessed$X)) {
      if (verbose) message("reduce batch_size temporarily")
      batch_end <- nrow(preprocessed$X)
      # reduced batch size
    } else {
      # regular batch_size
      batch_end <<- batch_start + batch_size  # 12 
    }
    if (verbose) message(paste("generating bach number",
                               batch_num, batch_start, "-", batch_end,
                               "of file", file, "index", next_file))
    x_batch <- preprocessed$X[batch_start:batch_end, , ]# dim shoiuld be (batch_size, length, words)
    y_batch <- preprocessed$Y[batch_start:batch_end,] #  # dim should be (batch_size, words)
    # return the file
    list(x_batch, y_batch)
  }
}
