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

#' Get cell states from a fasta file
#'
#' @param model keras model
#' @param fasta.path path to fasta file
#' @param maxlen time steps to unroll for
#' @param batch.size how many subsequences are predicted in parallel
#' @param verbose print output
#' @export
getStatesFromFasta <- function(model = NULL,
                               fasta.path = "example_files/fasta/a.fasta",
                      maxlen = 80,
                      batch.size = 100,
                      verbose = T){
  if (verbose)
    message("preprocess...")  
  # prepare fasta
  preprocessed <- deepG::preprocessFasta(fasta.path,
                           maxlen = maxlen,
                           vocabulary = c("\n", "a", "c", "g", "t"))
  batch.num <- 1
  batch.start <- 1
  batch.end <- batch.start + batch.size
  states <- NULL
  states <- list()
  while (batch.start < nrow(preprocessed$X)) {
    if ((batch.start + batch.size) > nrow(preprocessed$X)) {
      if (verbose)
        message("reduce batch.size temporarily")
      batch.end <- nrow(preprocessed$X)
      # reduced batch size
    }
    if (verbose)
      message(
        paste(
          "generating batch number",
          batch.num,
          batch.start,
          "-",
          batch.end
        ))
    x.batch <-
      preprocessed$X[batch.start:batch.end, , ] # dim shoiuld be (batch_size, length, words)
    states[[batch.num]] <- keras::predict_on_batch(model, x.batch)
    # update batch index
    batch.num <- batch.num + 1 
    batch.start <- batch.end + 1
    batch.end <- batch.start + batch.size
  }
  states.matrix <- do.call(rbind, states)
  return(states.matrix)
}
