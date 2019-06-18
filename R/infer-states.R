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
getStates <- function(model.path,
                      x,
                      maxlen = 30,
                      batch.size = 100,
                      run.name = "output",
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
  if (batch.size < 1)
    ArgumentCheck::addError(
      msg = "'batch.size' should be a positive integer",
      argcheck = Check
    )
  #* Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(Check)

  model <- load_model_hdf5(model.path)
  # Remove the last 2 layers
  keras::pop_layer(model)
  keras::pop_layer(model)
  states <- predict(model, x, batch_size = batch.size)
  # we dont have predictions in the beginning so create some empty cell response
  # so we set it to zero
  states.begining <- states[1:maxlen,] * 0
  states.final <- rbind(states.begining, states)
  # save states as hdf5
  if (verbose) print("saving states...")
  if (type == "hdf5") {
    file <- hdf5r::H5File$new(paste0(run.name, "_states.hdf5"), mode = "a")
    file.grp <- hdf5r::file.h5$create_group("states1")
    file.grp <- states.final
    hdf5r::h5close(file)
  } else {
    write.table(states.final,
                file = paste0(run.name, "_states.csv"),
                sep = ";", quote = F, col.names = F,
                row.names = F)

  }
  return(states.final)
}
