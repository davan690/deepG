#' Limit Tensorflow proccess to one or multiple GPUs.
#'
#' @param gpus GPU ids
#' @export
startGPUSession <- function(gpus = "0"){
  require(tensorflow)
  tf$reset_default_graph()
  sess_config <- list()
  Sys.setenv(CUDA_VISIBLE_DEVICES = gpus)
  sess_config$device_count <- list(GPU = 1L, CPU = 1L)
  session_conf <- do.call(tf$ConfigProto, sess_config)
  sess <- tf$Session(graph = tf$get_default_graph(), config = session_conf)
  sess$run(tf$global_variables_initializer())
  return(sess)
}

#' Ends the GPU sesstion
#'
#' @export
end_gpu_session <- function(){
  sess$close()
}

list.local.devices <- function(){
  message(tensorflow::tf$python$client$device_lib$list_local_devices())
}

gpu.available <- function() {
  check <- tensorflow::tf$python$test$is_gpu_available()
  if (check) {
    message("GPUs found!") 
  } else {
    message("No GPUs found!")
  }
  return(check)
}
