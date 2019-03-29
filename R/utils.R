#' generates hdf5 file containing the character id for each time step
#'
#' @param dat character
#' @param filename filename where hdf5 file is written to
#' @param verbose TRUE/FALSE
#' @export
writehdf5 <- function(dat, filename = "train.hdf5", verbose = F) {
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
writedict <- function(dat, filename = "train.dict") {
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