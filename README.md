![Coverage
Status](https://img.shields.io/codecov/c/github/hiddengenome/altum/master.svg)

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

```
data(crispr_sample) # load dataset named crispr_sample to workspace

# check if the vocabulary is ok, for DNA sequences only {A,C,G,T} and
# one new-genome character should be there
get_vocabulary(crispr_sample)

# reformat the training set as hdf5 char-ids
writehdf5(crispr_sample, file="train.hdf5")

# write dictionary
writedict(crispr_sample, file="train.dict")

# this generating one-hot representations of 30-character sized chunks and returns
# a list with elements X and Y storing training and target representations
crispr_sample_preprocessed <- preprocess(crispr_sample)

# train the model, will be saved as hdf5
history <- train_lstm(crispr_sample_preprocessed) # takes about 7 minutes on CPU

# extract states based on the model exported on the step before
getstates("example_run_full_model.hdf5", crispr_sample_preprocessed$X)
```