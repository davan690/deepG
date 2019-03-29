#' Get cell states of semi-redundant chunks
#'
#' @param model_path path to keras model in hdf5 format
#' @param x semi-redundant chunks (one-hot)
#' @param maxlen time steps to unroll for
#' @param batch_size how many samples are trained in parallel
#' @param run_name name of output files without ending
#' @param type will save as hdf5 if type is set to 'hdf5', otherwise as csv
#' @param verbose TRUE/FALSE
#' @export
getstates <- function(model_path,
                      x,
                      maxlen = 30,
                      batch_size = 100,
                      run_name = "output",
                      type = "csv",
                      verbose = F){
  require(dplyr)
  require(hdf5r)
  require(keras)

  Check <- ArgumentCheck::newArgCheck()
  #* Add an error if maxlen <1
  if (maxlen < 1)
    ArgumentCheck::addError(
      msg = "'maxlen' must be >= 1",
      argcheck = Check
    )
  #* Add an error if batch_size negative
  if (batch_size < 1)
    ArgumentCheck::addError(
      msg = "'batch_size' should be a positive integer",
      argcheck = Check
    )
  #* Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(Check)

  model <- load_model_hdf5(model_path)
  # Remove the last 2 layers
  keras::pop_layer(model)
  keras::pop_layer(model)
  states <- predict(model, x, batch_size = batch_size)
  # we dont have predictions in the beginning so create some empty cell response
  # so we set it to zero
  states_begining <- states[1:maxlen,] * 0
  final_states <- rbind(states_begining, states)
  # save states as hdf5
  if (verbose) print("saving states...")
  if (type == "hdf5") {
    file <- hdf5r::H5File$new(paste0(run_name, "_states.hdf5"), mode = "a")
    file.grp <- hdf5r::file.h5$create_group("states1")
    file.grp <- final_states
    hdf5r::h5close(file)
  } else {
    write.table(final_states,
                file = paste0(run_name, "_states.csv"),
                sep = ";", quote = F, col.names = F,
                row.names = F)

  }
  return(final_states)
}


#' Iterates over entries in NCBI list and export cell state responses from a model
#'
#' @param list_path path to file downloaded from NCBI
#' @export
getstates_ncbi <- function(list_path, model_path,  ...) {
  # process additional arguments
  dots = list(...)
  # load NCBI list
  ncbi_list <- read.table(list_path, header = F, sep = ",")
  for (i in 1:nrow(ncbi_list)) {
    preprocessed_genome <- altum::download_from_ncbi(ncbi_list[i,6])
    def.vals = list(model = model_path,
                    x = preprocessed_genome,
                    maxlen = 30,
                    run_name = paste0("batch_output_", ncbi_list[i,1]),
                    batch_size = 200,
                    type = "csv")
    ind = unlist(lapply(dots[names(def.vals)], is.null))
    dots[names(def.vals)[ind]] = def.vals[ind]
    states <- do.call(altum::getstates, dots)
  }
}