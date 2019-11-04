#' generates hdf5 file containing the character id for each time step
#'
#' @param dat character
#' @param filename filename where hdf5 file is written to
#' @param verbose TRUE/FALSE
#' @export
writeHdf5 <- function(dat, filename = "train.hdf5", verbose = F) {
	require(dplyr)
	require(plyr)
	require(hdf5r)
	tokenized <- dat %>%
		tokenizers::tokenize_characters(strip_non_alphanum = FALSE,
																		simplify = TRUE)
	# get sorted vocabulary
	charset <- sort(unique(unlist(dat_tokenized)))

	# get corresponding character numbers
	charset_index <- 1:length(charset)
	names(charset_index) <- charset

	# replace character by character index
	tokenized_index <- mapvalues(tokenized, from = charset,
															 to = charset_index)

	if (verbose) print("saving states...")
	file <- hdf5r::H5File$new(filename, mode = "a")
	file.grp <- hdf5r::file.h5$create_group("words")
	file.grp <- array(tokenized_index)
	hdf5r::h5close(file)
}


#' generates dictionary for LSTMVis
#'
#' @param dat character
#' @param dat filename where dict file is written to
#' @export
writeDict <- function(dat, filename = "train.dict") {
	require(dplyr)
	require(plyr)

	tokenized <- dat %>%
		tokenizers::tokenize_characters(strip_non_alphanum = FALSE,
																		simplify = TRUE)
	# get sorted vocabulary
	charset <- sort(unique(unlist(dat_tokenized)))

	dict <- data.frame(char = charset,
										 index = 1:length(charset))
	write.table(dict, file = filename, quote = F, row.names = F,
							col.names = F, sep = " ")
}

#' wrapper for message(sprintf)
#'
#' @export
messagef <- function (..., .newline = TRUE) 
{
  message(sprintf(...), appendLF = .newline)
}

#' prints TF version
#'
#' @export
print.tf.version <- function() {
  message(paste("Tensoflow", tensorflow::tf$`__version__`, "found."))
}

ensure.loaded <- function() {
  invisible(tensorflow::tf$`__version__`)
}

#' Calculate the steps per epoch
#' 
#' Do one full preprocessing iteration to the FASTA file to figure out what the observed
#' steps_per_epoch value is.
#' @param dir Input directory where .fasta files are located
#' @param batch.size Number of samples  
#' @param maxlen Length of the semi-redundant sequences
#' @param format File format
#' @param verbose TRUE/FALSE
#' @export
calculateStepsPerEpoch <-
  function(dir,
           batch.size = 256,
           maxlen = 250,
           format = "fasta",
           verbose = F) {
    library(xfun)
    library(Biostrings)
    steps.per.epoch <- 0
    fasta.files <- list.files(
      path = xfun::normalize_path(dir),
      pattern = paste0("*.", format),
      full.names = TRUE
    )
    for (file in fasta.files) {
      fasta.file <- Biostrings::readDNAStringSet(file)
      seq <- paste0(paste(fasta.file, collapse = "\n"), "\n")
      steps.per.epoch <-
        steps.per.epoch + ceiling((nchar(seq) - maxlen - 1) / batch.size)
    }
    return(steps.per.epoch)
}