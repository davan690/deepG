# deepG <img src="man/figures/logo.png" width="131px" height="140px" align="right" style="padding-left:10px;background-color:white;" />


[![Build Status](https://travis-ci.org/hiddengenome/deepG.svg?branch=master)](https://travis-ci.org/hiddengenome/deepG)
[![codecov](https://codecov.io/gh/hiddengenome/deepG/branch/master/graph/badge.svg)](https://codecov.io/gh/hiddengenome/deepG)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL%20v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)

## Overview

deepG is a package for generating LSTM models from genomic text and provides scripts for various common tasks such as the extraction of cell responses. It also comes with example datasets of genomic and human readable languages for testing.

## Installation & Usage

Please see our [Wiki](https://github.com/hiddengenome/deepG/wiki) for further installation instructions. It covers also usage instructions for multi-GPU machines
- [Installation on a desktop machine](https://github.com/hiddengenome/deepG/wiki/Installation-of-deepG-on-desktop)
- [Training of GenomeNet](https://github.com/hiddengenome/deepG/wiki/Howto-train-GenomeNet)

See the help files `?deepG` to get startet. 

## Datasets

The library comes with three different datasets for testing. 

- The set `data(parenthesis)` contains 1M characters of the parenthesis synthetic language generated from a very simple counting language with a parenthesis and letter alphabet Σ = {( ) 0 1 2 3 4 }. The language is constrained to match parentheses, and nesting is limited to at most 4 levels deep. Each opening parenthesis increases and each closing parenthesis decreases the nesting level, respectively. Numbers are generated randomly, but are constrained to indicate the nesting level at their position.  
- The set `data(crispr_full)` containging all CRISPR loci found in NCBI representative genomes with neighbor nucleotides up and downstream.
- The set `data(crispr_sample)` containging a subset of `crispr_full`

## License and copyright
Copyright 2019 Philipp Münch

Source code to deepG is made available under the terms of the GNU Affero General Public License (AGPL). deepG is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.

## Supported by

<p float="left">
  <img src="man/figures/hzi.jpg" width="200" />
  <img src="man/figures/dfg.jpg" width="200" />
  <img src="man/figures/bmbf.jpeg" width="100" /> 
  <img src="man/figures/aws.png" width="100" /> 
</p>
