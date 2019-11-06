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
  
  #if (verbose)
  #  print(sprintf("corpus length: %d", xt))
  
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
    dataset <- purrr::map(seq(1, length(text) - maxlen , by = 1),
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
      # generate one-hot encoding for one subset
      x[i, ,] <- sapply(vocabulary, function(x) {
        as.integer(x == dataset$sentece[[i]])
      })
      # target (next nucleotide in sequence)
      y[i,] <- as.integer(vocabulary == dataset$next_char[[i]])
    }
    } else {
      # use labels
      dataset <- purrr::map(seq(1, length(text) - maxlen , by = 1),
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
                            vocabulary = c("-", "|", "a", "c", "g", "t"),
                            verbose = F) {
  # process labels
  if (!is.null(labels)){
    fasta.file.labels <- Biostrings::readDNAStringSet(labels)
    seq.labels <- paste0("|", paste(fasta.file, collapse = "-"),"|") 
  } 
  
  # process corpus
  fasta.file <- Biostrings::readDNAStringSet(path)
  seq <- paste0("|", paste(fasta.file, collapse = "-"),"|") 
  
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


#' Helper function for fastaFileGenerator
#' @param sequence character sequence 
#' @param maxlen length of one sample
#' 
#' Returns one hot encoding for every sequence  
#' For example: sequence = "acatg", maxlen = 4, vocabulary = c("-", "a", "c", "g", "t") leads to
#' X = (0 1 0 0 0  
#'      0 0 1 0 0
#'      0 1 0 0 0
#'      0 0 0 0 1)
#' Y = (0 0 0 1 0)       
#' @export
sequenceToArray <- function(sequence, maxlen, vocabulary = c("-", "|", "a", "c", "g", "t")){
  len_voc <- length(vocabulary)
  len_seq <- nchar(sequence)
  z <- array(0L, dim=c(len_seq*len_voc))  
  tokenizer <- keras::fit_text_tokenizer(keras::text_tokenizer(char_level = TRUE, lower = TRUE), vocabulary) 
  sequence_int <- keras::texts_to_sequences(tokenizer, sequence) 
  seq_unlist <- sequence_int[[1]]
  adjust <- len_voc*(seq(0, len_seq - 1))
  
  # every row in z one-hot encodes one character in sequence
  z[adjust + seq_unlist] <- 1L
  z <- keras::array_reshape(z, dim=c(len_seq, len_voc))
  
  x <- array(0L, dim = c(len_seq - maxlen, maxlen, len_voc))
  for (i in 1:(len_seq - maxlen)){
    x[i, , ] <- z[i:(maxlen+i-1), ] 
  }
  y <- z[(maxlen+1):nrow(z),]
  list(x, y)
}


#' custom generator for fasta files, will produce chunks in size of batch.size
#' by iterating over the input files. 
#' @param corpus.dir input directory where .fasta files are located
#' @param format file format
#' @param batch.size number of samples  
#' @param maxlen length of one sample 
#' @param max_iter stop after max_iter number of iterations failed to produce new sample 
#' @param verbose TRUE/FALSE show information about files and running time  
#' @param seqStart insert character at beginning of sequence
#' @param seqEnd insert character at end of sequence
#' @param withinFile insert characters within sequence
#' @param randomFiles TRUE/FALSE, whether to go through files randomly or sequential 
#' @param showWarnings TRUE/FALSE, give warning if character outside vocabulary appears   
#' @export
fastaFileGenerator <- function(corpus.dir,
                               format = "fasta",
                               batch.size = 256,
                               maxlen = 250,
                               max_iter = 20,
                               seqStart = "|",
                               seqEnd= "|",
                               withinFile = "-",
                               vocabulary = c("|","-","a", "c", "g", "t"),
                               verbose = FALSE,
                               randomFiles = FALSE,
                               #replaceInFileSampling = TRUE,
                               #step = 1, 
                               showWarnings = TRUE){
  
  for (i in c(seqStart, seqStart, withinFile)) {
    if(!(i %in% vocabulary))
      stop("seqStart, seqEnd and withinFile variables must be in vocabulary")
  }  
  
  fasta.files <- list.files(
    path = xfun::normalize_path(corpus.dir),
    pattern = paste0("*.", format),
    full.names = TRUE)
  num_files <- length(fasta.files)
  
  # regular expression for chars outside vocabulary
  pattern <- paste0("[^", paste0(vocabulary, collapse = ""), "]")
  
  # sequence vector collects strings until one batch can be created 
  sequence_vector <- vector("character")
  sequence_vector_index <- 1
  
  if (randomFiles) {
    file_index <- sample(num_files, size = 1)
  } else {
    file_index <- 1
  }
  start_index <- 1
  next_sub_seq <- 1
  num_samples <- 0
  
  # pre-load the first file
  file <- fasta.files[[file_index]]
  fasta.file <- Biostrings::readDNAStringSet(file)
  seq <- paste0(seqStart, paste(fasta.file, collapse = withinFile), seqEnd)  
  # test for chars outside vocabulary
  if (showWarnings){
    charsOutsideVoc <- stringr::str_detect(seq, pattern)  
    if (charsOutsideVoc) warning("file ", file, " contains characters outside of vocabulary")
  }
  
  # split seq around chars not in vocabulary
  sub_seq <- splitSequence(seq = seq, vocabulary = vocabulary, maxlen = maxlen)
  current_seq <- sub_seq[1]
  length_current_seq <- nchar(current_seq)
  num_sub_seq <- length(sub_seq)
  
  if (verbose) message("initializing")
  
  function() {
    start_time <- Sys.time()
    iter <- 1
    # loop until enough samples collected
    while(num_samples < batch.size){  
      if (num_sub_seq > 1 & showWarnings) warning("file ", file, " contains characters outside of vocabulary")
      # loop through sub-sequences/files until sequence of suitable length is found   
      while(start_index + maxlen >= length_current_seq | is.na(length_current_seq)){
        next_sub_seq <<- next_sub_seq + 1
        start_index <<- 1
        
        # go to next file (if condition true) or next sub-sequence (else) 
        if (next_sub_seq > length(sub_seq)){
          next_sub_seq <<- 1
          if (randomFiles) {
            file_index <<- sample(num_files, size = 1)
          } else {
            file_index <<- file_index + 1
          }
          if (file_index > length(fasta.files)) file_index <<- 1
          file <<- fasta.files[[file_index]]
          fasta.file <<- Biostrings::readDNAStringSet(file)
          seq <- paste0(seqStart, paste(fasta.file, collapse = withinFile), seqEnd)  
          sub_seq <- splitSequence(seq = seq, vocabulary = vocabulary, maxlen = maxlen)
          current_seq <<- sub_seq[1]
          length_current_seq <<- nchar(current_seq) 
          num_sub_seq <<- length(sub_seq)
          # test for chars outside vocabulary
          if (showWarnings){
            charsOutsideVoc <- stringr::str_detect(seq, pattern)  
            if (charsOutsideVoc) warning("file ", file, " contains characters outside the vocabulary")
          }
        } else {
          current_seq <<- sub_seq[next_sub_seq]
          length_current_seq <<- nchar(current_seq) 
        } 
        
        if(iter > max_iter){
          stop('exceeded max_iter value, try reducing maxlen parameter')
          break
        }
        iter <- iter + 1
      }
      
      # go to end of sub sequence or stop when enough samples are collected 
      end_index <- min(start_index + maxlen + (batch.size - num_samples) - 1,
                       length_current_seq)
      current_sub_seq <- substr(current_seq, start_index, end_index)
      sequence_vector[sequence_vector_index] <- current_sub_seq
      sequence_vector_index <- sequence_vector_index + 1
      length_cur_sub_seq <- nchar(current_sub_seq)
      num_samples <- num_samples + length_cur_sub_seq - maxlen 
      start_index <<- start_index + length_cur_sub_seq - maxlen 
    }
    
    # one hot encode strings collected in sequence_vector and connect arrays
    array_list <- purrr::map(1:length(sequence_vector), ~sequenceToArray(sequence_vector[.x], maxlen = maxlen))
    x <- array_list[[1]][[1]]
    y <- array_list[[1]][[2]]
    if (length(array_list)>1){
      for (i in 2:length(array_list)){
        x <- abind::abind(x,array_list[[i]][[1]], along = 1)
        y <- abind::abind(y,array_list[[i]][[2]], along = 1)
      }
    }
    end_time <- Sys.time()
    if (verbose){
      cat("running time:", end_time - start_time, "\n") 
      obj_size <- format(object.size(list(x,y)),  units = "auto")
      cat("batch size:", obj_size, "\n")
    }  
    list(X = x, Y = y)
  }
}

