#' Example training dataset consiting of a sequence of nucleotides of CRISPR loci
#' Filtered for unambigous characters and contains only characters in the vocabulary {A,G,G,T
#' }
#' Can be loaded to workspace via data(crispr_sample)
#' @format Large character of 4.9 Mb
#' @references \url{http://github.com/philippmuench}
#' @examples
#' head(crispr_sample)
"crispr_sample"

#' Example training dataset consiting of a sequence of nucleotides of CRISPR loci
#' Filtered for unambigous characters and contains only characters in the vocabulary {A,G,G,T
#' }
#' contain all CRISPR loci found in NCBI representative genomes using CRT
#' Can be loaded to workspace via data(crispr_full)
#' @format Large character of 4.9 Mb
#' @references \url{http://github.com/philippmuench}
#' @examples
#' head(crispr_full)
"crispr_full"


#' Training dataset of synthetic parenthesis language
#' Can be loaded to workspace via data(parenthesis)
#' @format Large character of 4.9 Mb
#' @references \url{http://github.com/philippmuench}
#' @examples
#' head(parenthesis)
"parenthesis"

#' Get vocabulary from character string
#'
#' Use this function with a character string as input such as data(train)
#'
#' @param char a character string of text, length of one
#' @export
get_vocabulary <- function(char) {
  library(dplyr)
  vocabulary <-  char %>%  stringr::str_to_lower() %>%
    stringr::str_c(collapse = "\n") %>%
    tokenizers::tokenize_characters(strip_non_alphanum = FALSE, simplify = TRUE) %>%
    unique() %>%
    sort()
  print(vocabulary)
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
#' @export
#'
preprocess <- function(char, maxlen = 30, vocabulary) {
  require(dplyr)
  Check <- ArgumentCheck::newArgCheck()
  #* Add an error if maxlen <1
  if (maxlen < 1)
    ArgumentCheck::addError(
      msg = "'maxlen' must be >= 1",
      argcheck = Check
    )

  if(!missing(vocabulary)) {
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
  print(sprintf("corpus length: %d", length(text)))

  if(missing(vocabulary)) {
    vocabulary <- text %>%
      unique() %>%
      sort()
  }
  print(sprintf("vocabulary size: %d", length(vocabulary)))

  # Cut the text in semi-redundant sequences of maxlen characters
  print("generation of semi-redundant sequences ...")
  dataset <- purrr::map(seq(1, length(text) - maxlen - 1, by = 1),
                        ~ list(sentece = text[.x:(.x + maxlen - 1)],
                               next_char = text[.x + maxlen]))
  dataset <- purrr::transpose(dataset)

  # Vectorization
  x <-
    array(0, dim = c(length(dataset$sentece), maxlen, length(vocabulary)))
  y <-
    array(0, dim = c(length(dataset$sentece), length(vocabulary)))
  print("vectorization ...")
  pb <-
    txtProgressBar(min = 0,
                   max = length(dataset$sentece),
                   style = 3)
  for (i in 1:length(dataset$sentece)) {
    setTxtProgressBar(pb, i)
    x[i, ,] <- sapply(vocabulary, function(x) {
      as.integer(x == dataset$sentece[[i]])
    })
    y[i,] <- as.integer(vocabulary == dataset$next_char[[i]])
  }
  results <- list("X" = x, "Y" = y)
  return(results)
}

#' Train LSTM model
#'
#' @param x list with X and Y
#' @param run_name for saving purpose
#' @param maxlen time steps to unroll for (e.g. length of semi-redundant chunks)
#' @param dropout_rate dropout rate for LSTM
#' @param layer_size number of cells per network layer
#' @param batch_size how many chunks are trained in parallel
#' @param validation_split proportion of training set that will be used for validation
#' @param learning_rate learning rate for optimizer
#' @param layer_size number of cells per network layer#'
#' @param cudnn if true, using layer_cudnn_lstm() instead of layer_lstm()
#' @param vocabulary_size number of unique chars in training set'
#' @param epochs number of full iterations over the dataset
#' @export
train_lstm <- function(dat,
                        run_name = "run",
                        maxlen = 30,
                        dropout_rate = .3,
                        layer_size = 128,
                        batch_size = 128,
                        validation_split = 0.5,
                        learning_rate = 0.001,
                        cudnn = FALSE,
                        vocabulary_size = 5,
                        epochs = 1) {
  require(dplyr)
  require(h5)

  Check <- ArgumentCheck::newArgCheck()
  #* Add an error if maxlen <1
  if (maxlen < 1)
    ArgumentCheck::addError(
      msg = "'maxlen' must be >= 1",
      argcheck = Check
    )
  #* Add an error if dropout_rate is < 0 or > 1
  if (dropout_rate > 1 |  dropout_rate < 0)
    ArgumentCheck::addError(
      msg = "'dropout_rate' must be between 0 and 1",
      argcheck = Check
  )
  #* Add an error if layer_size negative
  if (layer_size < 1)
    ArgumentCheck::addError(
      msg = "'layer_size' should be a positive integer",
      argcheck = Check
  )
  #* Add an error if layer_size negative
  if (batch_size < 1)
    ArgumentCheck::addError(
      msg = "'batch_size' should be a positive integer",
      argcheck = Check
  )
  #* Add an error if dropout_rate is < 0 or > 1
  if (validation_split > 1 |  validation_split < 0)
    ArgumentCheck::addError(
      msg = "'validation_split' must be between 0 and 1",
      argcheck = Check
    )
  #* Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(Check)

  if (cudnn) {
    model <- keras::keras_model_sequential()
    model %>%
      keras::layer_cudnn_lstm(
        layer_size,
        input_shape = c(maxlen, vocabulary_size),
        return_sequences = T
      ) %>%
      keras::layer_dropout(rate = dropout_rate) %>%
      keras::layer_cudnn_lstm(layer_size) %>%
      keras::layer_dropout(rate = dropout_rate) %>%
      keras::layer_dense(vocabulary_size) %>%
      keras::layer_activation("softmax")
  } else {
    model <- keras::keras_model_sequential()
    model %>%
      keras::layer_lstm(
        layer_size,
        input_shape = c(maxlen, vocabulary_size),
        return_sequences = T
      ) %>%
      keras::layer_dropout(rate = dropout_rate) %>%
      keras::layer_lstm(layer_size) %>%
      keras::layer_dropout(rate = dropout_rate) %>%
      keras::layer_dense(vocabulary_size) %>%
      keras::layer_activation("softmax")
  }

  optimizer <- keras::optimizer_adam(lr = learning_rate)
  model %>% keras::compile(loss = "categorical_crossentropy",
                           optimizer = optimizer)

  history <- model %>% keras::fit(
    dat$X,
    dat$Y,
    batch_size = batch_size,
    validation_split = validation_split,
    epochs = epochs
  )

  optimizer <- keras::optimizer_adam(lr = learning_rate)
  # save training metadata
  print("save model...")
  Rmodel <- keras::serialize_model(model, include_optimizer = TRUE)
  save(Rmodel, file = paste0(run_name, "_full_model.Rdata"))
  keras::save_model_hdf5(
    model,
    paste0(run_name, "_full_model.hdf5"),
    overwrite = TRUE,
    include_optimizer = TRUE
  )
  return(history)
}

#' Get cell state responses of a model to a preprocessed text corpus
#'
#' @param model LSTM model in hdf5 format
#' @param x preprocessed dataset for prediction
#' @param maxlen time steps to unroll for (e.g. length of semi-redundant chunks)
#' @param batch_size how many chunks are trained in parallel
#' @param run_name name of output files without ending
#' @param type will save as hdf5 if type is set to 'hdf5', otherwise as csv
#' @export
getstates <- function(model_path,
                      x,
                      maxlen = 30,
                      batch_size = 100,
                      run_name = "output",
                      type = "csv"){
  require(dplyr)
  require(h5)
  require(keras)

  Check <- ArgumentCheck::newArgCheck()
  #* Add an error if maxlen <1
  if (maxlen < 1)
    ArgumentCheck::addError(
      msg = "'maxlen' must be >= 1",
      argcheck = Check
  )
  #* Add an error if batch_size negative
  if (batch_size < 1)
    ArgumentCheck::addError(
      msg = "'batch_size' should be a positive integer",
      argcheck = Check
  )
  #* Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(Check)

  model <- load_model_hdf5(model_path)
  # Remove the last 2 layers
  keras::pop_layer(model)
  keras::pop_layer(model)
  states <- predict(model, x, batch_size = batch_size)
  # we dont have predictions in the beginning so create some empty cell response
  # so we set it to zero
  states_begining <- states[1:maxlen,] * 0
  final_states <- rbind(states_begining, states)
  # save states as hdf5
  print("saving states...")
  if (type == "hdf5") {
    file <- h5::h5file(paste0(run_name, "_states.hdf5"), mode = "a")
    file["states1"] <- final_states
    h5::h5close(file)
  } else {
    write.table(final_states,
                file = paste0(run_name, "_states.csv"),
                sep = ";", quote = F, col.names = F,
                row.names = F)

  }
}