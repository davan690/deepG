  # extracts the start positions of all potential samples (considering step size and vocabulary)
  extractStartIndices <- function(fasta.path, maxlen, step, vocabulary, seqStart, withinFile, seqEnd){
    
    fasta.file <- Biostrings::readDNAStringSet(fasta.path)
    seq <- paste0(seqStart, paste(fasta.file, collapse = withinFile), seqEnd) %>% stringr::str_to_lower()
    
    len_seq <- nchar(seq)
    stopifnot(len_seq >= maxlen)
    # regular expressions for allowed characters
    voc_pattern <- paste0("[^", paste0(vocabulary, collapse = ""), "]")
    
    pos_of_amb_nucleotides <- stringr::str_locate_all(seq, pattern = voc_pattern)[[1]][,1]
    non_start_index <-  pos_of_amb_nucleotides - maxlen + 1
    
    # define range of unallowed start indices   
    non_start_index <- purrr::map(non_start_index, ~(.x:(.x + maxlen - 1))) %>% 
      unlist() %>% union((len_seq - maxlen + 1):len_seq) %>% unique() # TODO: go to end or where last prediction possible 
    # drop non-positive values 
    if (length(non_start_index[non_start_index < 1])){
      non_start_index <- unique(c(1, non_start_index[non_start_index >= 1]))
    }
    
    allowed_start <- setdiff(1:len_seq, non_start_index)
    len_start_vector <- length(allowed_start) 
    
    # only keep indices with sufficient distance, as defined by step 
    if (len_start_vector < 1) {
      message("Can not extract a single sampling point with current settings.")
      return(NULL)
    }
    index <- allowed_start[1]
    start_indices <- vector("integer")
    start_indices[1] <- allowed_start[1]
    count <- 1
    for (i in 1:len_start_vector-1){
      if (allowed_start[i + 1] - index >=step){
        count <- count + 1  
        start_indices[count] <- allowed_start[i + 1]
        index <- allowed_start[i + 1]
      }
    }
    start_indices
  }
  
  
  # extract samples (as strings) from one fasta file   
  sequenceList <- function(fasta.path, maxlen, step, vocabulary, seqStart, withinFile, seqEnd){
    fasta.file <- Biostrings::readDNAStringSet(fasta.path)
    seq <- paste0(seqStart, paste(fasta.file, collapse = withinFile), seqEnd) %>% stringr::str_to_lower()
    start_indices <- extractStartIndices(fasta.path = fasta.path, maxlen = maxlen, step = step, vocabulary = vocabulary,
                                         seqStart = seqStart, withinFile = withinFile, seqEnd = seqEnd)
    purrr::map(start_indices, ~substr(seq, start = .x, stop = .x + maxlen - 1))
  }
  
  # one-hot-encoding of list of sequences
  arrayList <- function(seqList, vocabulary){
    maxlen <- nchar(seqList[[1]])
    lapply(seqList, FUN = function(seq){
      len_voc <- length(vocabulary)
      len_seq <- nchar(seq)
      tokenizer <- keras::fit_text_tokenizer(keras::text_tokenizer(char_level = TRUE, lower = TRUE), vocabulary) 
      sequence_int <- keras::texts_to_sequences(tokenizer, seq)
      seq_unlist <- sequence_int[[1]]
      adjust <- len_voc * (seq(0, len_seq - 1))
      x <- array(0, dim = maxlen * len_voc)
      x[adjust + seq_unlist] <- 1
      x <- keras::array_reshape(x, dim=c(len_seq, len_voc))
      x
      }
    )
  }
  
  # 
  predictionGenerator <- function(fasta.path, batch.size, maxlen, vocabulary, step, seqStart, withinFile, seqEnd){
    start <- 1
    seqList <- sequenceList(fasta.path = fasta.path, maxlen = maxlen, step = step, vocabulary = vocabulary,
                            seqStart = seqStart, withinFile = withinFile, seqEnd = seqEnd)
    names(seqList) <- 1:length(seqList)
    x <- array(0, dim = c(batch.size, maxlen, length(vocabulary)))
    function(){
      if ((start + batch.size - 1) < length(seqList)){
        seqSubList <- seqList[start:(start + batch.size - 1)]
        arraySubList <- arrayList(seqSubList, vocabulary)
        for (i in 1:batch.size){
          x[i, , ] <- arraySubList[[i]]
        }
        start <<- start + batch.size
      # last batch may be smaller than batch.size   
      } else {
        x <- array(0, dim = c(length(seqList) - start + 1, maxlen, length(vocabulary)))
        seqSubList <- seqList[start:length(seqList)]
        arraySubList <- arrayList(seqSubList, vocabulary)
        for (i in 1:length(arraySubList)){
          x[i, , ] <- arraySubList[[i]]
        } 
        # print("Finished one iteration")
        start <<- 1
      }
      list(x)
      }
  }

writeStatesToH5 <- function(model.path, layer_depth, fasta.path, round_digits = 2, h5_filename,
                        step, vocabulary, seqStart, withinFile, seqEnd, batch.size, verbose = TRUE){
  
  stopifnot(batch.size > 0)
  stopifnot(layer_depth > 0)
  stopifnot(!file.exists(paste0(h5_filename,".h5")) & !file.exists(h5_filename))
  
  # load model and file
  model <- keras::load_model_hdf5(model.path)
  fasta.file <- Biostrings::readDNAStringSet(fasta.path)
  seq <- paste0(seqStart, paste(fasta.file, collapse = withinFile), seqEnd) %>% stringr::str_to_lower()
  
  # extract maxlen
  maxlen <- model$input$shape[1] 
  
  start_indices <- extractStartIndices(fasta.path = fasta.path, maxlen = maxlen, step = step, vocabulary = vocabulary,
                                       seqStart = seqStart, withinFile = withinFile, seqEnd = seqEnd)
  
  num_of_layers <- length(model$layers)
  stopifnot(num_of_layers >= layer_depth)
  if (verbose){
    model
    cat("Original model has", num_of_layers, "layers")
  }
  
  if (num_of_layers - layer_depth !=0){
    for (i in 1:(num_of_layers - layer_depth)){
      keras::pop_layer(model)
    }
  }
  
  # extract number of neurons in last layer
  if (length(model$output$shape$dims) == 3){
    layer.size <- model$output$shape[2]
  } else {
    layer.size <- model$output$shape[1]
  }
  
  model
  
  # create h5 file to store states
  if (!stringr::str_detect(h5_filename, pattern = "\\.h5$")){
    filename <- paste0(h5_filename,".h5")
  }
  h5_file <- hdf5r::H5File$new(filename, mode = "w") 
  h5_file[["states"]] <- array(0, dim = c(0, layer.size)) 
  writer <- h5_file[["states"]]
  
  # how often to call the predictionGenerator function
  numberOfGeneratorCalls <- ceiling(length(start_indices)/batch.size)
  
  # generator for input 
  gen <- predictionGenerator(fasta.path = fasta.path, batch.size = batch.size, maxlen = maxlen, vocabulary = vocabulary, 
                             step = step, seqStart = seqStart, withinFile = withinFile, seqEnd = seqEnd)
  
  row <- 1
  for (i in 1:numberOfGeneratorCalls){
    input <- gen()[[1]]   
    current_batch_size <- dim(input)[1]
    activations <- keras::predict_on_batch(object = model, x = input)
    activations <- as.array(activations)
    # some layers give predictions for every char in sequence, 
    # since return_sequences = TRUE, filter last prediction 
    if (length(dim(activations)) == 3){
      activations <- activations[ , maxlen, ]
    }
    writer[row:(row + current_batch_size - 1), ] <- round(activations, digits = round_digits)
    row <- row + current_batch_size
  }
  # save start indices of sequences  
  #hdf5r::h5attr(writer, "start_indices") <- as.character(start_indices)
  h5_file$close_all()
}

#readStatesFromH5 <- function(h5_path, ){
#  h5_files <- hdf5r::H5File$new(h5_path, mode = "r")
#  read_lines <- h5_files[["states"]]
#  
#  h5_file$close_all()
#}
