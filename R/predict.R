#' Predict next nucleotides in sequence output as S4 class
#'
#' @param sequence input sequence, length should be in sync with the model
#' if length exceeds input_shape of model then only the right side of the
#' sequence will be used
#' @param model trained lstm
#' @param vocabulary ordered vocabulary of input sequence
#' @export
predict_next_nucleotide <- function(sequence,
                                    model,
                                    vocabulary = c("\n", "a", "c", "g", "t")){
  Check <- ArgumentCheck::newArgCheck()
    if (missing(sequence))
    stop("Need to specify sequence")
  if (missing(model))
    stop("Need to specify model")
  Check <- ArgumentCheck::newArgCheck()
  #* Add an error if sequence length is too small
  if (nchar(sequence) < model$input_shape[2] )
    ArgumentCheck::addError(
      msg = paste0("'sequence' must be of length >", model$input_shape[2]),
      argcheck = Check
    )
  #* Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(Check)
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


#' Replaces specific nucleotides in a sequence one by one
#' @param sequence input sequence
#' @param model trained lstm
#' @param char character that will be replaced
#' @param vocabulary ordered vocabulary of input sequence
#' @export
replace_char <- function(sequence,
                         model,
                         char = "X",
                         vocabulary = c("\n", "a", "c", "g", "t")){
  require(stringr)
  Check <- ArgumentCheck::newArgCheck()
  if (missing(sequence))
    stop("Need to specify sequence")
  if (missing(model))
    stop("Need to specify model")
  Check <- ArgumentCheck::newArgCheck()
  #* Add an error if sequence length is too small
  if (nchar(sequence) < model$input_shape[2] )
    ArgumentCheck::addError(
      msg = paste0("'sequence' must be of length >", model$input_shape[2]),
      argcheck = Check
    )
  #* Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(Check)

  while (str_detect(sequence, char)) {
    # get the position
    next_position <- stringr::str_locate_all(pattern = 'X', sequence)[[1]][[1]]
    # seed text for model is the most-right chunk of text
    # with size of model$input_shape[[2]]
    seed <- substr(sequence,
                   next_position - model$input_shape[[2]] - 1,
                   next_position - 1)
    prediction <- predict_next_nucleotide(seed, model, vocabulary)
    sequence <- paste0(prediction@solution,
                       substr(sequence, next_position + 1,
                              nchar(sequence)))
  }
  return(sequence)
}

#' One beam
#' @param sequence input sequence
#' @param model trained lstm
#' @param vocabulary ordered vocabulary of input sequence
#' @param num_alternatives number of best solutions
#' @export
one_beam <- function(sequence, model, vocabulary, num_alternatives = 2){
  # getting num_alternative best indexes in vocabulary
  solutions <- list()
  confidences <- list()
  prediction <- predict_next_nucleotide(sequence, model, vocabulary)
  for (nth_best in 1:num_alternatives){
    solutions[[nth_best]] <- paste0(sequence, vocabulary[Rfast::nth(prediction@alternative_probabilty, nth_best, descending = T, index.return = T)])
    confidences[[nth_best]] <- Rfast::nth(prediction@alternative_probabilty, nth_best, descending = T)
  }
  return(list(solutions, confidences))
}