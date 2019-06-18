[![Build Status](https://travis-ci.org/hiddengenome/altum.svg?branch=master)](https://travis-ci.org/hiddengenome/altum)
[![codecov](https://codecov.io/gh/hiddengenome/altum/branch/master/graph/badge.svg)](https://codecov.io/gh/hiddengenome/altum)

## Overview

altum is a package for generating LSTM models from genomic text and provides scripts for various common tasks such as the extraction of cell responses. It also comes with example datasets of genomic and human readable languages for testing.

## installation

### dependencies

For model saving, this package requires the hdf5 library. You can install it from source 

``` bash
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.5/src/hdf5-1.10.5.tar.gz
tar xvzf hdf5-1.10.5.tar.gz
cd hdf5-1.10.5
./configure --prefix=/usr --enable-fortran --enable-cxx
make
make check
sudo make install
```

and then

``` r
install.packages("hdf5r", configure.args="--with-hdf5=/usr/bin/h5cc")
```

## this library 

``` r
install.packages("devtools")
devtools::install_github("hiddengenome/altum")
```

### enable GPU support

On default, Keras will be installed without GPU support. To support GPUs reinstall Keras via

``` r
install.packages("keras")
keras::install_keras(tensorflow = "gpu")
```

# usage

## model generation via the fasta generator

To build a language model based on a collection of fasta files located in a folder, please run 

``` r
history <- train_lstm_generator("input_dir")
```

and if you are working on a big machine (DGX-1)

``` r
history <- train_lstm_generator(path = "input_dir", cudnn = T, multiple_gpu = T, gpu_num = 1:8, run_name= "GenomeNet", epochs = 100, steps_per_epoch = 10000)
```

See the `?train_lstm_generator` for further information how to setup file names and network size.

## model generation using data held in the RAM

#### load a example dataset

Load dataset named crispr_sample to workspace using `data(crispr_sample)`. Then check if the vocabulary is as exprected, for DNA sequences only {A,C,G,T} and one new-genome character should be there `get_vocabulary(crispr_sample)`

#### preprocessing

Now generate one-hot representations of 30-nucleotide (here) sized chunks. This fucntion returns a list with elements X and Y storing training and target representations `crispr_sample_preprocessed <- preprocess(crispr_sample)`

#### training
Then we can train the model, will be saved as hdf5 `history <- train_lstm(crispr_sample_preprocessed)` which takes about 7 minutes on CPU. You can plot training and validation error via `plot(history)`. On default, the model will be written to `run_full_model.hdf5`. 

### batch prediction
After we trained the CrisprNet we want to use the model and predict some bacterial genomes. For this we can use the `getstates_ncbi()` which takes a list that can be downloaded fomr the NCBI genome browser directly. The hidden state responses to each genome will be saved to the disk and can then be used by uisum for visualization. For example the function call would be `getstates_ncbi("tmp/ncbi.csv", "run_full_model.hdf5")`

## Datasets

The library comes with three different datasets for testing. 

- The set `data(parenthesis)` contains 100k charakters of the parenthesis synthetic language generated from a very simple counting language with a parenthesis and letter alphabet Σ = {( ) 0 1 2 3 4 }. The language is constrained to match parentheses, and nesting is limited to at most 4 levels deep. Each opening parenthesis increases and each closing parenthesis decreases the nesting level, respectively. Numbers are generated randomly, but are constrained to indicate the nesting level at their position.  
- The set `data(crispr_full)` containging all CRISPR loci found in NCBI representative genomes with neighbor nucleotides up and downstream.
- The set `data(crispr_sample)` containging a subset of `crispr_full`

### LSTMVis support 

these function produce output that can be used for LSTMVis

```
# reformat the training set as hdf5 char-ids
writehdf5(crispr_sample, file="train.hdf5")

# write dictionary
writedict(crispr_sample, file="train.dict")

# extract states based on the model exported on the step before
getstates("example_run_full_model.hdf5", crispr_sample_preprocessed$X)
```
