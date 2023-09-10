
<!-- README.md is generated from README.Rmd. Please edit that file -->

# adestr <a href='https://github.com/jan-imbi/adestr'><img src='man/figures/sticker.png' align="right" height="170" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/jan-imbi/adestr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jan-imbi/adestr/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/jan-imbi/adestr/branch/master/graph/badge.svg?token=ORYWTYOZPT)](https://app.codecov.io/gh/jan-imbi/adestr?branch=master)
[![License](https://img.shields.io/badge/license-MIT-blue)](https://github.com/jan-imbi/adestr/blob/master/LICENSE.md)
<!-- badges: end -->

This package implements methods to evaluate the performance
characteristics of various point and interval estimators for adaptive
two-stage designs with prespecified sample-size recalculation rules.
Further, it allows for evaluation of these estimators on real datasets,
and it implements methods to calculate p-values.

Currently, it works for designs objects which were produced by the
R-package `adoptr`, which calculates optimal design parameters adaptive
two-stage designs.

## Installation

You can install the development version of adestr by typing

``` r
remotes::install_github("https://github.com/jan-imbi/adestr")
```

into your R console.

## Information for reviewers

The scripts to reproduce the results from the paper can be found in the
`/data/code/` directory of this repository. The results themselves are
located in the `/data/` directory.

The easiest way to inspect the results is to [clone this
repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository).

## General example for usage of the package

Here is a quick example showing the capabilities of `adestr`. First,
load `adestr`:

``` r
library(adestr)
#> Loading required package: adoptr
```

Then, you can evaluate the performance of an estimator like this:

``` r
evaluate_estimator(
 score = MSE(),
 estimator = SampleMean(),
 data_distribution = Normal(two_armed = TRUE),
 design = get_example_design(),
 mu = c(0, 0.3, 0.6),
 sigma = 1
)
#> Design:                              TwoStageDesign<n1=28;0.8<=x1<=2.3;n2=10-40>
#> Data Distribution:                                             Normal<two-armed>
#> Estimator:                                                           Sample mean
#> Assumed sigma:                                                                 1
#> Assumed mu:                                                          0.0 0.3 0.6
#> Results:
#>  Expectation:                                -0.03523827  0.28169661  0.63556747
#>  Bias:                                       -0.03523827 -0.01830339  0.03556747
#>  Variance:                                      0.05558910 0.07330464 0.06591361
#>  MSE:                                           0.05683084 0.07363966 0.06717865
```

You can analyze a dataset like this:

``` r
set.seed(321)
dat <- data.frame(
 endpoint = c(rnorm(28, .2, 1), rnorm(28, 0, 1),
              rnorm(23, .2, 1), rnorm(23, 0, 1)),
 group = factor(rep(c("ctl", "trt", "ctl", "trt"),
                    c(28,28,23,23))),
 stage = rep(c(1L, 2L), c(56, 46))
)
analyze(
 data = dat,
 estimator = get_example_statistics(),
 data_distribution = Normal(),
 sigma = 1,
 design = get_example_design()
)
#> Design:                              TwoStageDesign<n1=28;0.8<=x1<=2.3;n2=10-40>
#> Data Distribution:                                             Normal<two-armed>
#> Z1                                                                          1.75
#> Actual number of stages:                                                       2
#> 
#> Stage 2 results:
#>  Sample mean:                                                          0.5044685
#>  Pseudo Rao-Blackwellized:                                             0.3559511
#>  Median unbiased (LR test ordering):                                   0.4806717
#>  Bias reduced MLE (iterations=1):                                      0.4994132
#>  SWCF ordering CI:                                       [0.06262736, 0.7232758]
#>  LR test ordering CI:                                     [0.2127897, 0.7949281]
#>  SWCF ordering p-value:                                                 0.010977
#>  LR test ordering p-value:                                          0.0001877094
```
