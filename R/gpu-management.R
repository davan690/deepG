#' Limit Tensorflow proccess to one or multiple GPUs.
#'
#' @param gpus number of GPU to use (as integer), where first GPU on machine is 1
#' @param growth memory consumption of GPU grows and will not be blocked
#' @param env set environmental variable CUDA_VISIBLE_DEVICES
#' @export
startGPUSession <- function(gpu.number, growth = T, env = F) {
  # please note that CUDA_VISIBLE_DEVICES prevents processes from doing cudaMemcpy from/to devices not owned by the process. There's a significant performance degradation when NCCL is used with P2P communication disabled.
  if (env)
    Sys.setenv(CUDA_VISIBLE_DEVICES = gpu.number) 

  # get available GPUs on the machine as a list
  available.gpus <- tensorflow::tf$config$experimental$get_visible_devices('GPU')
  if (gpu.number >= length(available.gpus)) {
    tensorflow::tf$config$experimental$set_visible_devices(available.gpus[[gpu.number]], 'GPU')
  } else {
    stop("gpu.number is out of bounds with GPUs on machine.")
  }
  if (growth)
    tensorflow::tf$config$experimental$set_memory_growth(device = gpu, enable = TRUE)
}

list.local.devices <- function(){
  message(tensorflow::tensorflow::tf$config$experimental$get_visible_devices('GPU'))
}

is.gpu.available <- function() {
  res <- tryCatch({
    tensorflow::tf$test$is_gpu_available()
  }, error = function(e) {
    warning("Can not determine if GPU is configured.", call. = FALSE);
    NA
  })
  res
}

is.cuda.build <- function() {
  res <- tryCatch({
    tensorflow::tf$test$is_built_with_cuda()
  }, error = function(e) {
    warning("Can not determine if TF is build with CUDA", call. = FALSE);
    NA
  })
  res
}
