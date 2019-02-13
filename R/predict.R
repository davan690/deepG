#' Predict next nucleotides in sequence output as S4 class
#'
#' @param sequence input sequence
#' @param model trained lstm
#' @param vocabulary ordered vocabulary of input sequence
#' @export
predict_next_nucleotide <- function(sequence,
                                    model,
                                    vocabulary = c("\n", "a", "c", "g", "t")){
  if (missing(sequence))
    stop("Need to specify sequence")
  if (missing(model))
    stop("Need to specify model")
  Check <- ArgumentCheck::newArgCheck()
  #* Add an error if sequence length is too small
  if (nchar(sequence) < 1)
    ArgumentCheck::addError(
      msg = "'sequence' must be of length > 0",
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
             alternative_probabilty = "matrix"
           )
  )
  sentence <- sequence %>% stringr::str_to_lower()  %>%
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
             alternative_probabilty = preds))
}