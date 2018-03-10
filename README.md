# phylomoments

[![Travis-CI Build Status](https://travis-ci.org/dunleavy005/phylomoments.svg?branch=master)](https://travis-ci.org/dunleavy005/phylomoments)

This R package implements a simulation-free dynamic programming algorithm for
calculating higher-order moments of stochastic mapping summaries on a
phylogeny [3]. Utility functions related to simulation-based stochastic mapping [1]
and first-order moments [2] are also provided.

## Installation

First, install the `devtools` package (if it isn't already installed) using `install.packages("devtools")`.
Check to make sure you have a working development environment by running `devtools::has_devel()`.
If it returns `TRUE`, then your development environment is correctly set up; otherwise, you'll need to install additional tools.
More information on this subject can be found at https://github.com/hadley/devtools.
In addition, if you're using Mac OS X, you will need to run the following commands in the Terminal.

```shell
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```

The `phylomoments` package can then be installed from GitHub using the `install_github` function.

```r
library(devtools)
install_github("dunleavy005/phylomoments")
```

By default, vignettes are not included in the installation.
To install the `phylomoments` package with vignettes, use `install_github("dunleavy005/phylomoments", build_vignettes = TRUE)`.
Note that installing `phylomoments` with vignettes is a little slower than installing the package without them.

## Vignettes

We provide three vignettes that illustrate how to use the functions provided in this package.

1. Introductory Examples
2. Rate Variation Testing
3. Evolutionary Conservation Testing

You can access the `phylomoments` vignettes with `browseVignettes("phylomoments")`.

## References

1. Nielsen R (2002) "Mapping mutations on phylogenies", *Systematic Biology*, 51(5):729-739.

2. Minin VN and Suchard MA (2008) "Counting labeled transitions in continuous-time Markov models of evolution", *Journal of Mathematical Biology*, 56(3):391-412.

3. Dhar A and Minin VN (2017) "Calculating Higher-Order Moments of Phylogenetic Stochastic Mapping Summaries in Linear Time", *Journal of Computational Biology*, 24(5):377-399.
