#' deepG for GenomeNet
#'
#' deepG is a is an open source software librar for building deep neuronal 
#' networks for genomic modeling
#' 
#' This package generates \href{http://www.genomenet.de}{GenomeNet}
#'
#' For additional documentation on the deepG package see
#' \href{https://genomenet.de}{https://genomenet.de}
#'
#' @import reticulate
#' @import keras
#' @import tensorflow
#'
#' @docType package
#' @name deepG
NULL

# globals
.globals <- new.env(parent = emptyenv())
.globals$tensorboard <- NULL

.onLoad <- function(libname, pkgname) {
  message("The deepG package has been successfully. Please see ?deepG for informations.")
}
