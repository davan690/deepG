#' Predict next nucleotides in sequence 
#' 
#' The output is a S4 class.
#'
#' @param sequence input sequence, length should be in sync with the model.
#' If length exceeds input.shape of model then only the right side of the
#' sequence will be used.
#' @param model trained model from the function \code{trainNetwork()}
#' @param vocabulary vocabulary of input sequence
#' @param verbose TRUE/FALSE
#' @examples 
#' \dontrun{
#' example.model <- keras::load_model_hdf5("example_model.hdf5")
#' sequence <- strrep("A", 100)
#' perdictNextNucleotide(sequence, example.model)}
#' @export
predictNextNucleotide <- function(sequence,
                                    model,
                                    vocabulary =  c("l", "a", "c", "g", "t"),
                                    verbose = F){
  
  stopifnot(!missing(sequence))
  stopifnot(!missing(model))
  stopifnot(nchar(sequence) >= model$input_shape[2])
  
  substringright <- function(x, n){
    substr(x, nchar(x)- n + 1, nchar(x))
  }
  # sequence can be longer then model input shape
  # if so just use the last input_shape chars
  
  sentence <- tokenizers::tokenize_characters(
    stringr::str_to_lower(substringright(sequence, as.numeric(model$input_shape[2]))),
    strip_non_alphanum = FALSE, simplify = TRUE)
  
  x <- sapply(vocabulary, function(x){
    as.numeric(x == sentence)
  })
  x <- keras::array_reshape(x, c(1, dim(x)))

  if(verbose) {
    message("Prediction ...")}
  
  preds <- keras::predict_proba(model, x)
  next_index <- which.max(preds)
  next_char <- vocabulary[next_index]
  # return a S4 class
  return(new("prediction",
             next_char = next_char,
             probability = preds[next_index],
             index = next_index,
             alternative_probability = preds,
             solution = paste0(sequence, next_char)))
}


#' Replaces specific nucleotides in a sequence
#' 
#' @param sequence input sequence, length should be in sync with the model.
#' If length exceeds input.shape of model then only the right side of the
#' sequence will be used.
#' @param model trained model from the function \code{trainNetwork()}
#' @param char character in the sequence that will be replaced
#' @param vocabulary ordered vocabulary of input sequence
#' @examples 
#' \dontrun{
#' example.model <- keras::load_model_hdf5("example_model.hdf5")
#' replaceChar(sequence = sequence, model = example.model)}
#' @export
replaceChar <- function(sequence,
                         model,
                         char = "X",
                         vocabulary =  c("l", "a", "c", "g", "t")){
  
  stopifnot(!missing(sequence))
  stopifnot(!missing(model))
  stopifnot(nchar(sequence) >= model$input_shape[2])

  while (stringr::str_detect(sequence, char)) {
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
