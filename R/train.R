#' Train LSTM model via fasta generator
#'
#' @param path path to fasta files
#' @param run_name name of the run (without file ending)
#' @param maxlen time steps to unroll for (e.g. length of semi-redundant chunks)
#' @param dropout_rate dropout rate for LSTM
#' @param layer_size number of cells per network layer
#' @param batch_size how many chunks are trained in parallel
#' @param num_layers_lstm number of LSTM layers
#' @param codon_cnn first layer is a CNN layer with size of 3 
#' @param validation_split proportion of training set that will be used for validation
#' @param learning_rate learning rate for optimizer
#' @param layer_size number of cells per network layer#'
#' @param cudnn if true, using layer_cudnn_lstm() instead of layer_lstm()
#' @param multiple_gpu if true, multi_gpu_model will be used
#' @param gpu_num number of GPUs to be used, only relevant if multiple_gpu is true
#' @param cpu_merge true on default, false recommend for NVlink, only relevant if multiple_gpu is true
#' @param vocabulary_size number of unique chars in training set'
#' @param epochs number of full iterations over the dataset
#' @param max_queue_size queue on fit_generator()
#' @param steps_per_epoch Total number of steps (batches of samples) to yield from generator before declaring one epoch finished and starting the next epoch. It should typically be equal to the number of samples if your dataset divided by the batch size
#' @param verbose TRUE/FALSE
#' @export
train_lstm_generator <- function(path,
                       run_name = "run",
                       maxlen = 30,
                       dropout_rate = .3,
                       layer_size = 128,
                       batch_size = 256,
                       num_layers_lstm = 2,
                       codon_cnn = T,
                       validation_split = 0.05,
                       learning_rate = 0.001,
                       cudnn = FALSE,
                       multiple_gpu = FALSE,
                       cpu_merge = TRUE,
                       gpu_num = 2,
                       vocabulary_size = 5,
                       epochs = 1,
                       max_queue_size = 50,
                       steps_per_epoch = 20,
                       verbose = F) {
  require(dplyr)
  library(magrittr)
  require(tensorflow)
  
  Check <- ArgumentCheck::newArgCheck()
  #* Add an error if maxlen <1
  if (maxlen < 1)
    ArgumentCheck::addError(msg = "'maxlen' must be >= 1",
                            argcheck = Check)
  #* Add an error if dropout_rate is < 0 or > 1
  if (dropout_rate > 1 |  dropout_rate < 0)
    ArgumentCheck::addError(msg = "'dropout_rate' must be between 0 and 1",
                            argcheck = Check)
  #* Add an error if layer_size is smaller than 1
  if (layer_size < 1)
    ArgumentCheck::addError(msg = "'layer_size' should be a positive integer",
                            argcheck = Check)
  #* Add an error if number of layers is smaller than 1
  if (num_layers_lstm < 1)
    ArgumentCheck::addError(msg = "'num_layers_lstm' should be a positive integer",
                            argcheck = Check)
  #* Add an error if layer_size negative
  if (batch_size < 1)
    ArgumentCheck::addError(msg = "'batch_size' should be a positive integer",
                            argcheck = Check)
  #* Add an error if dropout_rate is < 0 or > 1
  if (validation_split > 1 |  validation_split < 0)
    ArgumentCheck::addError(msg = "'validation_split' must be between 0 and 1",
                            argcheck = Check)
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
  
  if (cudnn) {
    for (i in 1:(num_layers_lstm - 1)) {
      model %>%
        keras::layer_cudnn_lstm(
          layer_size,
          input_shape = c(maxlen, vocabulary_size),
          return_sequences = T
        ) %>%
        keras::layer_dropout(rate = dropout_rate)
    }
    # last LSTM layer should be with return_sequences = F
    model %>%
      keras::layer_cudnn_lstm(layer_size) %>%
      keras::layer_dropout(rate = dropout_rate)
    
  } else {
    if (codon_cnn) {
      model %<>% 
        keras::layer_conv_1d(kernel_size = 3, # 3 aa are a codon
                             padding = "same",
                             activation = "relu",
                             filters = 81,
                             input_shape = c(maxlen, vocabulary_size))  %>%
        layer_max_pooling_1d(pool_size = 3)  %>%
        layer_batch_normalization(momentum = .8)
    }
    for (i in 1:(num_layers_lstm - 1)) {
      model %<>% 
        keras::layer_lstm(
          layer_size,
          input_shape = c(maxlen, vocabulary_size),
          return_sequences = T
        ) %>%
        keras::layer_dropout(rate = dropout_rate)
    }
    # last LSTM layer should be with return_sequences = F
    model %>%
      keras::layer_lstm(layer_size) %>%
      keras::layer_dropout(rate = dropout_rate)
  }
  
  model %>% keras::layer_dense(vocabulary_size) %>%
    keras::layer_activation("softmax")
  optimizer <- keras::optimizer_rmsprop(lr = learning_rate)
  
  if (multiple_gpu) {
    parallel_model <- keras::multi_gpu_model(model,
                                             gpus = gpu_num,
                                             cpu_merge = cpu_merge)
    parallel_model %>% keras::compile(loss = "categorical_crossentropy",
                                      optimizer = optimizer)
    gen <- fasta_files_generator(path)
    
    start.time <- Sys.time()
    
    # fit using generator
    history <- model %>% fit_generator(
      generator = gen,
      steps_per_epoch = 100, # will auto-reset after see all sample
      max_queue_size = 50,
      epochs = epochs
    )
    
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
    # set-up the fasta generator
    gen <- fasta_files_generator(path)
    start.time <- Sys.time()

    # fit using generator
    history <- model %>% fit_generator(
      generator = gen,
      steps_per_epoch = 100, # will auto-reset after see all sample
      max_queue_size = 50,
      epochs = epochs
    )
    end.time <- Sys.time()
  }
  
  print(paste("Training time:", end.time - start.time))
  
  # save training metadata
  if (verbose)
    print("save model...")
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


#' Train  LSTM model
#'
#' @param dat preprocessed input data (semi-redundant chunks)
#' @param run_name name of the run (without file ending)
#' @param maxlen time steps to unroll for (e.g. length of semi-redundant chunks)
#' @param dropout_rate dropout rate for LSTM
#' @param layer_size number of cells per network layer
#' @param batch_size how many chunks are trained in parallel
#' @param num_layers_lstm number of LSTM layers
#' @param codon_cnn first layer is a CNN layer with size of 3 
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
                       num_layers_lstm = 2,
                       codon_cnn = T,
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
  library(magrittr)
  require(tensorflow)
  
  Check <- ArgumentCheck::newArgCheck()
  #* Add an error if maxlen <1
  if (maxlen < 1)
    ArgumentCheck::addError(msg = "'maxlen' must be >= 1",
                            argcheck = Check)
  #* Add an error if dropout_rate is < 0 or > 1
  if (dropout_rate > 1 |  dropout_rate < 0)
    ArgumentCheck::addError(msg = "'dropout_rate' must be between 0 and 1",
                            argcheck = Check)
  #* Add an error if layer_size is smaller than 1
  if (layer_size < 1)
    ArgumentCheck::addError(msg = "'layer_size' should be a positive integer",
                            argcheck = Check)
  #* Add an error if number of layers is smaller than 1
  if (num_layers_lstm < 1)
    ArgumentCheck::addError(msg = "'num_layers_lstm' should be a positive integer",
                            argcheck = Check)
  #* Add an error if layer_size negative
  if (batch_size < 1)
    ArgumentCheck::addError(msg = "'batch_size' should be a positive integer",
                            argcheck = Check)
  #* Add an error if dropout_rate is < 0 or > 1
  if (validation_split > 1 |  validation_split < 0)
    ArgumentCheck::addError(msg = "'validation_split' must be between 0 and 1",
                            argcheck = Check)
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
  
  if (cudnn) {
    for (i in 1:(num_layers_lstm - 1)) {
      model %>%
        keras::layer_cudnn_lstm(
          layer_size,
          input_shape = c(maxlen, vocabulary_size),
          return_sequences = T
        ) %>%
        keras::layer_dropout(rate = dropout_rate)
    }
    # last LSTM layer should be with return_sequences = F
    model %>%
      keras::layer_cudnn_lstm(layer_size) %>%
      keras::layer_dropout(rate = dropout_rate)
    
  } else {
    if (codon_cnn) {
      model %<>% 
        keras::layer_conv_1d(kernel_size = 3, # 3 aa are a codon
                             padding = "same",
                             activation = "relu",
                             filters = 81,
                             input_shape = c(maxlen, vocabulary_size))  %>%
        layer_max_pooling_1d(pool_size = 3)  %>%
        layer_batch_normalization(momentum = .8)
    }
      for (i in 1:(num_layers_lstm - 1)) {
        model %<>% 
          keras::layer_lstm(
            layer_size,
            input_shape = c(maxlen, vocabulary_size),
            return_sequences = T
          ) %>%
          keras::layer_dropout(rate = dropout_rate)
      }
    # last LSTM layer should be with return_sequences = F
    model %>%
      keras::layer_lstm(layer_size) %>%
      keras::layer_dropout(rate = dropout_rate)
  }
  
  model %>% keras::layer_dense(vocabulary_size) %>%
    keras::layer_activation("softmax")
  optimizer <- keras::optimizer_rmsprop(lr = learning_rate)
  
  if (multiple_gpu) {
    parallel_model <- keras::multi_gpu_model(model,
                                             gpus = gpu_num,
                                             cpu_merge = cpu_merge)
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
  if (verbose)
    print("save model...")
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

