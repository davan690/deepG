library(readr)
library(h5)

example_model <- load_model_hdf5("example_model.hdf5")
use_data(example_model, overwrite = TRUE)