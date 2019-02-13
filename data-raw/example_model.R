library(readr)
library(keras)

example_model <- keras::load_model_hdf5("data-raw/example_model.hdf5")
use_data(example_model, overwrite = TRUE)