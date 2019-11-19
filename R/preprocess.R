#' Returns the vocabulary from character string
#'
#' Use this function with a character string.
#'
#' @param char Character string of text with the length of one
#' @param verbose TRUE/FALSE
#' @export
getVocabulary <- function(char, verbose = F) {
  
  stopifnot(!is.null(char))
  stopifnot(nchar(char) > 0)
  
  vocabulary <- sort(unique(tokenizers::tokenize_characters(
    stringr::str_c(stringr::str_to_lower(char),collapse = "\n"), strip_non_alphanum = FALSE, simplify = TRUE)))

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
#' @param char Character input string of text with the length of one
#' @param maxlen Length of the semi-redundant sequences
#' @param vocabulary Char contains the vocabulary from the input char, it should be sorted.
#' If no vocabulary exists, it is generated from the input char
#' @param verbose TRUE/FALSE
#' @example preprocessSemiRedudant("abcd",maxlen=2,verbose=F)
#' @export
preprocessSemiRedundant <- function(char,
                                    maxlen = 250,
                                    vocabulary,
                                    verbose = F) {
  
  stopifnot(!is.null(char))
  stopifnot(nchar(char) > 0)
  stopifnot(maxlen > 0)
  
  # Load, collapse, and tokenize text ("ACGT" -> "a" "c" "g" "t")
  text <- tokenizers::tokenize_characters(stringr::str_c(stringr::str_to_lower(char), collapse = "\n"), strip_non_alphanum = FALSE, simplify = TRUE)
  
  # Generating vocabulary from input char with the function getVocabulary()
  if (missing(vocabulary)) {
    if (verbose)
      message("Finding the vocabulary ...")
    vocabulary <- getVocabulary(char)
  }
  
  if(verbose)
    message("Vocabulary size:", length(vocabulary))
  # Cut the text in semi-redundant sequences of maxlen characters
  
  if (verbose)
    message("Generation of semi-redundant sequences ...")
  
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
  
  results <- list("X" = x, "Y" = y)
  return(results)
}

#' Wrapper of the preprocessSemiRedundant()-function 
#' 
#' It called on the genomic contents of one
#' FASTA file. Multiple entries are combined with newline characters.
#' @param path Path to the FASTA file
#' @param maxlen Length of the semi-redundant sequences
#' @param vocabulary Char contains the vocabulary from the input char, it should be sorted.
#' If no vocabulary exists, it is generated from the input char
#' @param verbose TRUE/FALSE
#' @export
preprocessFasta <- function(path,
                            maxlen = 250,
                            vocabulary = c("l", "p", "a", "c", "g", "t"),
                            verbose = F) {

  
  # process corpus
  fasta.file <- Biostrings::readDNAStringSet(path)
  seq <- paste0("l", paste(fasta.file, collapse = "p"),"l") 
  
  if(verbose)
    message("Preprocessing the data ...")
  
  seq.processed <-
    preprocessSemiRedundant(char = seq, maxlen = maxlen, vocabulary = vocabulary,
                            verbose = F) 
  return(seq.processed)
}


#' Helper function for fastaFileGenerator
#' @param sequence character sequence 
#' @param maxlen length of one sample
#' @param vocabulary set of characters to encode  
#' @param step how often to take a sample
#' Returns one hot encoding for every sequence  
#' For example: sequence = "acatg", maxlen = 4, vocabulary = c("p","a", "c", "g", "t"), step = 1  leads to
#' X = (0 1 0 0 0  
#'      0 0 1 0 0
#'      0 1 0 0 0
#'      0 0 0 0 1)
#' Y = (0 0 0 1 0)       
#' @export
sequenceToArray <- function(sequence, maxlen, vocabulary, step){
  len_voc <- length(vocabulary)
  len_seq <- nchar(sequence)
  # len_seq should be n * step + maxlen + 1 for some integer n
  # how many samples can be extracted from sequence
  numberOfSamples <- ((len_seq - maxlen - 1)/step) + 1  
  z <- array(0L, dim=c(len_seq * len_voc))  
  tokenizer <- keras::fit_text_tokenizer(keras::text_tokenizer(char_level = TRUE, lower = TRUE), vocabulary) 
  sequence_int <- keras::texts_to_sequences(tokenizer, sequence) 
  seq_unlist <- sequence_int[[1]]
  adjust <- len_voc*(seq(0, len_seq - 1))
  
  # every row in z one-hot encodes one character in sequence
  z[adjust + seq_unlist] <- 1L
  z <- keras::array_reshape(z, dim=c(len_seq, len_voc))
  
  x <- array(0L, dim = c(numberOfSamples, maxlen, len_voc))
  for (i in 1:numberOfSamples){
    start <- 1 + (i - 1) * step
    end <- start + maxlen - 1
    x[i, , ] <- z[start:end, ] 
  }
  y <- z[1 + step * (0:(numberOfSamples-1)) + maxlen, ]
  list(x, y)
}


#' Helper function for fastaFileGenerator 
#' 
#' splits a character sequences into vector of sequences
#' @param seq a character sequence
#' @param vocabulary vector of allowed characters
#' @param maxlen sequences in the output < (maxlen + 1) get discarded    
#' @export 
splitSequence <- function(seq, vocabulary, maxlen){
  seq <- stringr::str_to_lower(seq)
  voc_pattern = paste0("[^", paste0(vocabulary, collapse = ""), "]")
  subSeq  <- stringr::str_split(seq, voc_pattern)[[1]]
  # only keep sequences that are long enough for one sample 
  subSeq <- subSeq[nchar(subSeq) > maxlen] 
  subSeq
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
#' @param step how often to take a sample#'
#' @param showWarnings TRUE/FALSE, give warning if character outside vocabulary appears   
#' @export
fastaFileGenerator <- function(corpus.dir,
                               format = "fasta",
                               batch.size = 256,
                               maxlen = 250,
                               max_iter = 20,
                               seqStart = "l",
                               seqEnd= "l",
                               withinFile = "p",
                               vocabulary = c("l","p","a", "c", "g", "t"),
                               verbose = FALSE,
                               randomFiles = FALSE,
                               step = 1, 
                               showWarnings = FALSE){
  
  for (i in c(seqStart, seqStart, withinFile)) {
    if(!(i %in% vocabulary) & i!="")
      stop("seqStart, seqEnd and withinFile variables must be in vocabulary")
  }  
  
  fasta.files <- list.files(
    path = xfun::normalize_path(corpus.dir),
    pattern = paste0("*.", format),
    full.names = TRUE)
  num_files <- length(fasta.files)
  
  if (randomFiles) fasta.files <- sample(fasta.files, replace = FALSE)
  
  # regular expression for chars outside vocabulary
  pattern <- paste0("[^", paste0(vocabulary, collapse = ""), "]")
  
  # sequence vector collects strings until one batch can be created   
  sequence_vector <- vector("character")
  sequence_vector_index <- 1
  
  file_index <- 1
  start_index <- 1
  next_sub_seq <- 1
  num_samples <- 0
  
  # pre-load the first file
  filePath <- fasta.files[[file_index]]
  fasta.file <- Biostrings::readDNAStringSet(filePath)
  seq <- paste0(seqStart, paste(fasta.file, collapse = withinFile), seqEnd)  
  
  # split seq around chars not in vocabulary
  seq_split <- splitSequence(seq = seq, vocabulary = vocabulary, maxlen = maxlen)
  current_seq <- seq_split[1]
  length_current_seq <- nchar(current_seq)
  num_sub_seq <- length(seq_split)
  # test for chars outside vocabulary
  if (showWarnings){
    charsOutsideVoc <- stringr::str_detect(stringr::str_to_lower(seq), pattern)  
    if (charsOutsideVoc) warning("file ", filePath, " contains characters outside vocabulary")
  }
  
  if (verbose) message("initializing")
  
  function() {
    iter <- 1
    # loop until enough samples collected
    while(num_samples < batch.size){  
      # loop through sub-sequences/files until sequence of suitable length is found   
      while((start_index + maxlen > length_current_seq) | is.na(length_current_seq)){
        next_sub_seq <<- next_sub_seq + 1
        start_index <<- 1
        
        # go to next file (if condition true) or next sub-sequence (else) 
        if (next_sub_seq > length(seq_split)){
          next_sub_seq <<- 1
          file_index <<- file_index + 1 
          if (file_index > length(fasta.files)) file_index <<- 1
          filePath <<- fasta.files[[file_index]]
          fasta.file <<- Biostrings::readDNAStringSet(filePath)
          seq <<- paste0(seqStart, paste(fasta.file, collapse = withinFile), seqEnd)  
          seq_split <<- splitSequence(seq = seq, vocabulary = vocabulary, maxlen = maxlen)
          current_seq <<- seq_split[1]
          length_current_seq <<- nchar(current_seq) 
          num_sub_seq <<- length(seq_split)
          # test for chars outside vocabulary
          if (showWarnings){
            charsOutsideVoc <- stringr::str_detect(stringr::str_to_lower(seq), pattern)  
            if (charsOutsideVoc) warning("file ", filePath, " contains characters outside vocabulary")
          }
        } else {
          current_seq <<- seq_split[next_sub_seq]
          length_current_seq <<- nchar(current_seq) 
        } 
        if(iter > max_iter){
          stop('exceeded max_iter value, try reducing maxlen parameter')
          break
        }
        iter <- iter + 1
      }
      
      # go as far as possible in sub-sequence or stop when enough samples are collected 
      potential_num_samples <- ceiling((length_current_seq - start_index - maxlen + 1)/step)
      end_index <- min(start_index + maxlen + (batch.size - num_samples - 1) * step,
                       (potential_num_samples - 1)*step + start_index + maxlen)
      
      sample_sub_seq <- substr(current_seq, start_index, end_index)
      sequence_vector[sequence_vector_index] <- sample_sub_seq
      length_sample_sub_seq <- nchar(sample_sub_seq)
      num_new_samples <- ((length_sample_sub_seq - maxlen - 1 ) / step) + 1  
      num_samples <- num_samples + num_new_samples 
      start_index <<- start_index + step * num_new_samples  
      sequence_vector_index <- sequence_vector_index + 1
     }
    
    # one hot encode strings collected in sequence_vector and connect arrays
    array_list <- purrr::map(1:length(sequence_vector),
                             ~sequenceToArray(sequence_vector[.x], maxlen = maxlen, vocabulary = vocabulary,  step = step))
    x <- array_list[[1]][[1]]
    y <- array_list[[1]][[2]]
    if (length(array_list) > 1){
      for (i in 2:length(array_list)){
        x <- abind::abind(x, array_list[[i]][[1]], along = 1)
        y <- rbind(y, array_list[[i]][[2]])
      }
    }
    
    # empty sequence_vector for next batch 
    sequence_vector <<- vector("character")
    sequence_vector_index <<- 1
    num_samples <<- 0
    
    list(X = x, Y = y)
  }
}
