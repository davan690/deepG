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

#' Trained LSTM model (only for testing)
#' @format keras model
#' @references \url{http://github.com/philippmuench}
#' @examples
#' head(example_model)
"example_model"