[![Build Status](https://travis-ci.org/hiddengenome/altum.svg?branch=master)](https://travis-ci.org/hiddengenome/altum) ![Coverage
Status](https://img.shields.io/codecov/c/github/hiddengenome/altum.svg?token=e2914fef-6e65-4145-b931-03d31277099f)

## Overview

altum is a package for generating LSTM models from genomic text and provides scripts for various common tasks such as the extraction of cell responses. It also comes with example datasets of genomic and human readable languages for testing.

## installation

``` r
#check for remotes
if(!require(remotes))
  install.packages("remotes")

# install altum
remotes::install_github("hiddengenome/altum")
```

### GPU support

this packaes requires the R interface for keras. 

``` r
install.packages("keras")
library(keras)
install_keras(tensorflow = "gpu")
```

for computation on CPU you can use `install_keras()`

## Datasets

The library comes with three different datasets for testing. 

- The set `data(parenthesis)` contains 100k charakters of the parenthesis synthetic language generated from a very simple counting language with a parenthesis and letter alphabet Σ = {( ) 0 1 2 3 4 }. The language is constrained to match parentheses, and nesting is limited to at most 4 levels deep. Each opening parenthesis increases and each closing parenthesis decreases the nesting level, respectively. Numbers are generated randomly, but are constrained to indicate the nesting level at their position.  
- The set `data(crispr_full)` containging all CRISPR loci found in NCBI representative genomes with neighbor nucleotides up and downstream.
- The set `data(crispr_sample)` containging a subset of `crispr_full` for testing

## Usage

### generate CrisprNet

#### load dataset
load dataset named crispr_sample to workspace using `data(crispr_sample)`. Then check if the vocabulary is as exprected, for DNA sequences only {A,C,G,T} and one new-genome character should be there `get_vocabulary(crispr_sample)`

#### preprocessing

Now generate one-hot representations of 30-nucleotide (here) sized chunks. This fucntion returns a list with elements X and Y storing training and target representations `crispr_sample_preprocessed <- preprocess(crispr_sample)`

#### training
Then we can train the model, will be saved as hdf5 `history <- train_lstm(crispr_sample_preprocessed)` which takes about 7 minutes on CPU. You can plot training and validation error via `plot(history)`. On default, the model will be written to `run_full_model.hdf5`. 

### batch prediction
After we trained the CrisprNet we want to use the model and predict some bacterial genomes. For this we can use the `getstates_ncbi()` which takes a list that can be downloaded fomr the NCBI genome browser directly. The hidden state responses to each genome will be saved to the disk and can then be used by uisum for visualization. For example the function call would be `getstates_ncbi("tmp/ncbi.csv", "run_full_model.hdf5")`


### training full CrisprNet
To train the full _CrisprNet_ model we use `data(crispr_full)` and following prarameter configuration in `train_lstm()` on a V100 GPU.

| Parameter        | Value  |
| ---------------- | ------ |
| maxlen           | 50     |
| layer_size       | 512    |
| batch_size       | 512    |
| learning_rate    | 0.0005 |
| validation_split | 0.05   |
| cudnn            | `TRUE` |
| epochs           | 50     |

Code to train the CrisprNet model:

```
library(altum)
data(crispr_full)
crispr_preprocessed <- preprocess(crispr_full, maxlen = 50) # needs up to 25G RAM
start_gpu_session() # limit process to one GPU
history <- train_lstm(crispr_preprocessed, run_name = "CrisprNet", maxlen = 50, layer_size = 512, batch_size = 512, learning_rate = 0.0005, cudnn = TRUE, epochs = 1, validation_split = 0.05) # needs up to 60G RAM
save(history, file="CrisprNet_history.Rdata")
end_gpu_session()
```

Code to get cell states from genomes
```
library(altum)
start_gpu_session(gpus="1") # limit process to one GPU
getstates_ncbi("ncbi.csv", "CrisprNet_full_model_test.hdf5")
```


### LSTMVis

these function produce output that can be used for LSTMVis

```
# reformat the training set as hdf5 char-ids
writehdf5(crispr_sample, file="train.hdf5")

# write dictionary
writedict(crispr_sample, file="train.dict")

# extract states based on the model exported on the step before
getstates("example_run_full_model.hdf5", crispr_sample_preprocessed$X)
```
