#' Predict next nucleotides in sequence 
#' 
#' The output is a S4 class.
#'
#' @param sequence Input sequence, length should be in sync with the model.
#' If length exceeds input.shape of model then only the right side of the
#' sequence will be used.
#' @param model Trained model from the function \code{trainNetwork()}
#' @param vocabulary Ordered vocabulary of input sequence
#' @export
predictNextNucleotide <- function(sequence,
                                    model,
                                    vocabulary = c("\n", "a", "c", "g", "t")){
  
  stopifnot(!missing(sequence))
  stopifnot(!missing(model))
  stopifnot(nchar(sequence) >= model$input_shape[2])
  
  require(keras)
  require(dplyr)
  require(tokenizers)
  
  setClass(Class = "prediction",
           representation(
             next_char = "character",
             probability = "numeric",
             index = "numeric",
             alternative_probabilty = "matrix",
             solution = "character"
           )
  )
  substringright <- function(x, n){
    substr(x, nchar(x)- n + 1, nchar(x))
  }
  # sequence can be longer then model input shape
  # if so just use the last input_shape chars
  sentence <- substringright(sequence, as.numeric(model$input_shape[2])) %>%
    stringr::str_to_lower() %>%
    tokenizers::tokenize_characters(strip_non_alphanum = FALSE, simplify = TRUE)
  x <- sapply(vocabulary, function(x){
    as.integer(x == sentence)
  })
  x <- array_reshape(x, c(1, dim(x)))

  preds <- keras::predict_proba(model, x)
  next_index <- which.max(preds)
  next_char <- vocabulary[next_index]
  # retrun a S4 class
  return(new("prediction",
             next_char = next_char,
             probability = preds[next_index],
             index = next_index,
             alternative_probabilty = preds,
             solution = paste0(sequence, next_char)))
}


#' Replaces specific nucleotides in a sequence
#' 
#' @param sequence Input sequence, length should be in sync with the model.
#' If length exceeds input.shape of model then only the right side of the
#' sequence will be used.
#' @param model Trained model from the function \code{trainNetwork()}
#' @param char Character in the sequence that will be replaced
#' @param vocabulary Ordered vocabulary of input sequence
#' @export
replaceChar <- function(sequence,
                         model,
                         char = "X",
                         vocabulary = c("\n", "a", "c", "g", "t")){
  require(stringr)
  
  stopifnot(!missing(sequence))
  stopifnot(!missing(model))
  stopifnot(nchar(sequence) >= model$input_shape[2])

  while (str_detect(sequence, char)) {
    # get the position
    next_position <- stringr::str_locate_all(pattern = 'X', sequence)[[1]][[1]]
    # seed text for model is the most-right chunk of text
    # with size of model$input_shape[[2]]
    seed <- substr(sequence,
                   next_position - model$input_shape[[2]] - 1,
                   next_position - 1)
    prediction <- predictNextNucleotide(seed, model, vocabulary)
    sequence <- paste0(prediction@solution,
                       substr(sequence, next_position + 1,
                              nchar(sequence)))
  }
  return(sequence)
}
