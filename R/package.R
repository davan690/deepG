#' deepG
#'
#' DeepG is an open source software library for generating and processing
#' \href{http://www.genomenet.de}{GenomeNet} a deep neuronal network for 
#' genomic modeling. This library allows you to train GenomeNet on a multiple 
#' GPU machine
#''
#' For additional documentation on the deepG package see
#' \href{http://www.genomenet.de}{http://www.genomenet.de}
#'
#' @import keras
#' @import tensorflow
#' @import reticulate
#' @import dplyr
#'
#' @docType package
#' @name deepG
NULL

# globals
.globals <- new.env(parent = emptyenv())
.globals$tensorboard <- NULL

.onLoad <- function(libname, pkgname) {
  # ensure that tensorflow initializes
  ensure_loaded()
  message("DeepG loaded!")
}
