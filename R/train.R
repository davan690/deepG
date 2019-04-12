#' Train standart LSTM model "Training time: 5.21544417142868" (2gpus)
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
#' @param multiple_gpu if true, multi_gpu_model will be used
#' @param gpu_num number of GPUs to be used, only relevant if multiple_gpu is true
#' @param cpu_merge true on default, false recommend for NVlink, only relevant if multiple_gpu is true
#' @param vocabulary_size number of unique chars in training set'
#' @param epochs number of full iterations over the dataset#
#' @param verbose TRUE/FALSE
#' @export
train_lstm <- function(dat,
                       run_name = "run",
                       maxlen = 30,
                       dropout_rate = .3,
                       layer_size = 128,
                       batch_size = 256,
                       validation_split = 0.05,
                       learning_rate = 0.001,
                       cudnn = FALSE,
                       multiple_gpu = FALSE,
                       cpu_merge = TRUE,
                       gpu_num = 2,
                       vocabulary_size = 5,
                       epochs = 1,
                       verbose = F) {
  require(dplyr)
  require(tensorflow)

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
  n#* Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(Check)

  # initialize model
  if (multiple_gpu) {
    # init template model under a CPU device scope
    with(tf$device("/cpu:0"), {
      model <- keras::keras_model_sequential()
    })
  } else {
    model <- keras::keras_model_sequential()
  }

  if (cudnn ) {
    model %>%
      keras::layer_cudnn_lstm(
        layer_size,
        input_shape = c(maxlen, vocabulary_size),
        return_sequences = T
      ) %>%
      keras::layer_dropout(rate = dropout_rate) %>%
      keras::layer_cudnn_lstm(layer_size) %>%
      keras::layer_dropout(rate = dropout_rate)

  } else {
    model %>%
      keras::layer_lstm(
        layer_size,
        input_shape = c(maxlen, vocabulary_size),
        return_sequences = T
      ) %>%
      keras::layer_dropout(rate = dropout_rate) %>%
      keras::layer_lstm(layer_size) %>%
      keras::layer_dropout(rate = dropout_rate)
  }

  model %>% keras::layer_dense(vocabulary_size) %>%
    keras::layer_activation("softmax")
  optimizer <- keras::optimizer_rmsprop(lr = learning_rate)

  if (multiple_gpu) {
    parallel_model <- keras::multi_gpu_model(model,
                                             gpus = gpu_num,
                                             cpu_merge= cpu_merge)
    parallel_model %>% keras::compile(loss = "categorical_crossentropy",
                             optimizer = optimizer)
    start.time <- Sys.time()
    history <- parallel_model %>% keras::fit(
      dat$X,
      dat$Y,
      batch_size = batch_size,
      validation_split = validation_split,
      epochs = epochs
    )
    end.time <- Sys.time()
  } else {
    model %>% keras::compile(loss = "categorical_crossentropy",
                             optimizer = optimizer)
    start.time <- Sys.time()
    history <- model %>% keras::fit(
      dat$X,
      dat$Y,
      batch_size = batch_size,
      validation_split = validation_split,
      epochs = epochs
    )
    end.time <- Sys.time()
  }

  print(paste("Training time:", end.time - start.time))

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