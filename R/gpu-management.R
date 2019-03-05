#' Set # of GPUs to be used in TensorFlow backend.
#'
#' @param gpus GPU ids
#' @export
#'
start_gpu_session <- function(gpus = "0"){
  require(tensorflow)
  tf$reset_default_graph()
  sess_config <- list()
  Sys.setenv(CUDA_VISIBLE_DEVICES = gpus)
  sess_config$device_count <- list(GPU = 1L, CPU = 1L)
  session_conf <- do.call(tf$ConfigProto, sess_config)
  sess <- tf$Session(graph = tf$get_default_graph(), config = session_conf)
  sess$run(tf$global_variables_initializer())
}

end_gpu_session <- function(){
  sess$close()
}
