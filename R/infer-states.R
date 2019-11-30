#' Get cell states of semi-redundant chunks
#'
#' @param model.path path to keras model in hdf5 format
#' @param x semi-redundant chunks (one-hot)
#' @param maxlen time steps to unroll for
#' @param batch.size how many samples are trained in parallel
#' @param run.name name of output files without ending
#' @param type save-type, will save as hdf5 if type is set to 'hdf5' (default .csv)
#' @param verbose TRUE/FALSE
#' @export

getStates <- function(model.path,
                      x,
                      maxlen = 30,
                      batch.size = 100,
                      run.name = "output",
                      type = "csv",
                      verbose = F){
  
  stopifnot(maxlen > 0)
  stopifnot(batch.size > 0)

  model <- keras::load_model_hdf5(model.path)
  # Remove the last 2 layers
  keras::pop_layer(model)
  keras::pop_layer(model)
  states <- keras::predict(model, x, batch_size = batch.size)
  # we dont have predictions in the beginning so create some empty cell response (set it zero)
  states.begining <- states[1:maxlen,] * 0
  states.final <- rbind(states.begining, states)
  # save states as hdf5
  if (verbose) print("Saving states ...")
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
                               verbose = F){
  if (verbose)
    message("Preprocess ...")  
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
        message("Reduce batch.size temporarily")
      batch.end <- nrow(preprocessed$X)
      # reduced batch size
    }
    if (verbose)
      message(
        paste(
          "Generating batch number",
          batch.num,
          batch.start,
          "-",
          batch.end
        ))
    x.batch <-
      preprocessed$X[batch.start:batch.end, , ] # dim should be (batch_size, length, words)
    states[[batch.num]] <- keras::predict_on_batch(model, x.batch)
    # update batch index
    batch.num <- batch.num + 1 
    batch.start <- batch.end + 1
    batch.end <- batch.start + batch.size
  }
  states.matrix <- do.call(rbind, states)
  return(states.matrix)
}
