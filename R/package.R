#' deepG for GenomeNet
#'
#' deepG is a is an open source software library for building deep neuronal 
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
#' @import dplyr
#'
#' @docType package
#' @name deepG
NULL

# globals
.globals <- new.env(parent = emptyenv())
.globals$tensorboard <- NULL

.onLoad <- function(libname, pkgname) {
  
  ensure.loaded <- function(x) {
    invisible(tensorflow::tf$`__version__`)
  }
  
  packageStartupMessage("The deepG package has been successfully loaded. Please see ?deepG for informations or the Wiki to get started https://github.com/hiddengenome/deepG/wiki")

  # check for GPU support
  if (is.gpu.available() & is.cuda.build()) {
    packageStartupMessage("To use GPUs, please run startGPUSession()")
  } else if (is.gpu.available() & !is.cuda.build()) {
    packageStartupMessage("GPUs are available, but Tensorflow is not supporting CUDA, please see Wiki (ttps://github.com/hiddengenome/deepG/wiki) for installation instructions")
    } else {
      packageStartupMessage("GPUs not found - deepG runs without GPU support!")
    }
  }
