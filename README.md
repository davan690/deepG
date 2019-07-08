[![Build Status](https://travis-ci.org/hiddengenome/deepG.svg?branch=master)](https://travis-ci.org/hiddengenome/deepG)
[![codecov](https://codecov.io/gh/hiddengenome/deepG/branch/master/graph/badge.svg)](https://codecov.io/gh/hiddengenome/deepG)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL%20v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)

## Overview

deepG is a package for generating LSTM models from genomic text and provides scripts for various common tasks such as the extraction of cell responses. It also comes with example datasets of genomic and human readable languages for testing.

## Installation

### Dependencies

For model saving, this package requires the `hdf5r` library which is based on hdf5 binaries. You can install it from source using the following commands. 

``` bash
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.5/src/hdf5-1.10.5.tar.gz
tar xvzf hdf5-1.10.5.tar.gz
cd hdf5-1.10.5
./configure --prefix=/usr --enable-fortran --enable-cxx
make
make check
sudo make install
```

``` r
install.packages("hdf5r", configure.args="--with-hdf5=/usr/bin/h5cc")
```

## Library 

``` r
install.packages("devtools")
devtools::install_github("hiddengenome/deepG")
```

### Enable GPU support

On default, Keras will be installed without GPU support. To support GPUs reinstall Keras via

``` r
install.packages("keras")
keras::install_keras(tensorflow = "gpu")
```

on AWS instances the GPU version can be installed via `conda create -y --name r-tensorflow tensorflow-gpu python=3.6.8` 

# Usage

## generating of GenomeNet 

### Generate a genomic language model from a collection of FASTA files

``` r
history <- trainNetwork(path = "input_dir")
```

When you have multiple GPUs available you can run

``` r
history <- trainNetwork(path = "input_dir", cudnn = T, multiple_gpu = T, gpu_num = 1:8, run_name= "GenomeNet", epochs = 100, steps_per_epoch = 10000)
```

See the `?trainNetwork` for further information how to setup file names and network size.

### Generate a genomic language model using data held in the RAM

``` r
# load example dataset
data(crispr_sample)

# generate one-hot encoding of 30-nucleotide sized chunks
crispr_preprocessed <- preprocessSemiRedundant(crispr_sample)
history <- trainNetwork(dataset = crispr_preprocessed)
```


## Train a deep network for multiclass prediction

here the FASTA files are located as individual files in `input_dir/` and labels are defined in `fasta-like` files in `fasta_label_dir/` with a `.txt` ending.

example file `file1.fasta` in `input_dir/`:

```
>fasta_header
ACGTGAGAGAGAGACAGAGATAGACAGAGATTATAA
```

example label file `file1.txt` in `fasta_label_dir/`, where {A,B,C} are the class labels.

```
>fasta_header
AAAAAAABBBBBBBAAAAAACCCCCCCCCCCCCCCC
```

``` r
history <- trainNetwork(path = "input_dir", labels = "fasta_label_dir/", label.vocabulary.size = 3)
```

## Datasets

The library comes with three different datasets for testing. 

- The set `data(parenthesis)` contains 100k charakters of the parenthesis synthetic language generated from a very simple counting language with a parenthesis and letter alphabet Σ = {( ) 0 1 2 3 4 }. The language is constrained to match parentheses, and nesting is limited to at most 4 levels deep. Each opening parenthesis increases and each closing parenthesis decreases the nesting level, respectively. Numbers are generated randomly, but are constrained to indicate the nesting level at their position.  
- The set `data(crispr_full)` containging all CRISPR loci found in NCBI representative genomes with neighbor nucleotides up and downstream.
- The set `data(crispr_sample)` containging a subset of `crispr_full`

## License and copyright
Copyright 2019 Philipp Münch

Source code to deepG is made available under the terms of the [GNU Affero General Public License](LICENSE.txt) (AGPL). deepG is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.

## Supported by

<p float="left">
  <img src="man/figures/hzi.jpg" width="200" />
  <img src="man/figures/dfg.jpg" width="200" />
  <img src="man/figures/bmbf.jpeg" width="100" /> 
  <img src="man/figures/aws.png" width="100" /> 
</p>