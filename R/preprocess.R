#' Get vocabulary from character string
#'
#' Use this function with a character string as input such as data(train)
#'
#' @param char a character string of text, length of one
#' @export
get_vocabulary <- function(char) {
  library(dplyr)
  stopifnot(!is.null(char))
  stopifnot(nchar(char)>0)
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