% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyze.R
\name{analyze}
\alias{analyze}
\alias{analyze,data.frame-method}
\title{Analyze a dataset}
\usage{
analyze(
  data,
  estimator,
  data_distribution,
  use_full_twoarm_sampling_distribution = FALSE,
  design,
  sigma,
  exact = FALSE
)

\S4method{analyze}{data.frame}(
  data,
  estimator,
  data_distribution,
  use_full_twoarm_sampling_distribution = FALSE,
  design,
  sigma,
  exact = FALSE
)
}
\arguments{
\item{data}{a data.frame containing the data to be analyzed}

\item{estimator}{object of class \code{PointEstimator}, \code{IntervalEstimator} or \code{PValue}.}

\item{data_distribution}{object of class \code{Normal} or \code{Student}.}

\item{use_full_twoarm_sampling_distribution}{logical indicating whether this estimator is intended to be used
with the full sampling distribution in a two-armed trial}

\item{design}{object of class \code{TwoStageDesign}.}

\item{sigma}{assumed standard deviation.}

\item{exact}{logical indicating usage of exact n2 function.}
}
\value{
\code{Results} object containing the values of the estimators
when applied to data.
}
\description{
Analyze a dataset
}
\examples{
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
  estimator = c(get_example_estimators()),
  data_distribution = Normal(),
  sigma = 1,
  design = get_example_design()
)
}