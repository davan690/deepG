#' Train standart LSTM model
#'
#' @param dat preprocessed input data (semi-redundant chunks)
#' @param run_name name of the run (without file ending)
#' @param maxlen time steps to unroll for (e.g. length of semi-redundant chunks)
#' @param dropout_rate dropout rate for LSTM
#' @param layer_size number of cells per network layer
#' @param batch_size how many chunks are trained in parallel
#' @param validation_split proportion of training set that will be used for validation
#' @param learning_rate learning rate for optimizer
#' @param layer_size number of cells per network layer#'
#' @param cudnn if true, using layer_cudnn_lstm() instead of layer_lstm()
#' @param vocabulary_size number of unique chars in training set'
#' @param epochs number of full iterations over the dataset#
#' @param verbose TRUE/FALSE
#' @export
train_lstm <- function(dat,
                        run_name = "run",
                        maxlen = 30,
                        dropout_rate = .3,
                        layer_size = 128,
                        batch_size = 128,
                        validation_split = 0.05,
                        learning_rate = 0.001,
                        cudnn = FALSE,
                        vocabulary_size = 5,
                        epochs = 1,
                       verbose = F) {
  require(dplyr)
  require(h5)

  Check <- ArgumentCheck::newArgCheck()
  #* Add an error if maxlen <1
  if (maxlen < 1)
    ArgumentCheck::addError(
      msg = "'maxlen' must be >= 1",
      argcheck = Check
    )
  #* Add an error if dropout_rate is < 0 or > 1
  if (dropout_rate > 1 |  dropout_rate < 0)
    ArgumentCheck::addError(
      msg = "'dropout_rate' must be between 0 and 1",
      argcheck = Check
  )
  #* Add an error if layer_size negative
  if (layer_size < 1)
    ArgumentCheck::addError(
      msg = "'layer_size' should be a positive integer",
      argcheck = Check
  )
  #* Add an error if layer_size negative
  if (batch_size < 1)
    ArgumentCheck::addError(
      msg = "'batch_size' should be a positive integer",
      argcheck = Check
  )
  #* Add an error if dropout_rate is < 0 or > 1
  if (validation_split > 1 |  validation_split < 0)
    ArgumentCheck::addError(
      msg = "'validation_split' must be between 0 and 1",
      argcheck = Check
    )
  #* Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(Check)

  if (cudnn) {
    model <- keras::keras_model_sequential()
    model %>%
      keras::layer_cudnn_lstm(
        layer_size,
        input_shape = c(maxlen, vocabulary_size),
        return_sequences = T
      ) %>%
      keras::layer_dropout(rate = dropout_rate) %>%
      keras::layer_cudnn_lstm(layer_size) %>%
      keras::layer_dropout(rate = dropout_rate) %>%
      keras::layer_dense(vocabulary_size) %>%
      keras::layer_activation("softmax")
  } else {
    model <- keras::keras_model_sequential()
    model %>%
      keras::layer_lstm(
        layer_size,
        input_shape = c(maxlen, vocabulary_size),
        return_sequences = T
      ) %>%
      keras::layer_dropout(rate = dropout_rate) %>%
      keras::layer_lstm(layer_size) %>%
      keras::layer_dropout(rate = dropout_rate) %>%
      keras::layer_dense(vocabulary_size) %>%
      keras::layer_activation("softmax")
  }

  optimizer <- keras::optimizer_adam(lr = learning_rate)
  model %>% keras::compile(loss = "categorical_crossentropy",
                           optimizer = optimizer)

  history <- model %>% keras::fit(
    dat$X,
    dat$Y,
    batch_size = batch_size,
    validation_split = validation_split,
    epochs = epochs
  )

  optimizer <- keras::optimizer_adam(lr = learning_rate)
  # save training metadata
  if (verbose) print("save model...")
  Rmodel <- keras::serialize_model(model, include_optimizer = TRUE)
  save(Rmodel, file = paste0(run_name, "_full_model.Rdata"))
  keras::save_model_hdf5(
    model,
    paste0(run_name, "_full_model.hdf5"),
    overwrite = TRUE,
    include_optimizer = TRUE
  )
  return(history)
}