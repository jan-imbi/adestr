% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimators.R, R/neyman_pearson.R
\docType{class}
\name{IntervalEstimator-class}
\alias{IntervalEstimator-class}
\alias{IntervalEstimator}
\alias{VirtualIntervalEstimator}
\alias{RepeatedCI}
\alias{StagewiseCombinationFunctionOrderingCI}
\alias{MLEOrderingCI}
\alias{LikelihoodRatioOrderingCI}
\alias{ScoreTestOrderingCI}
\alias{NaiveCI}
\alias{NeymanPearsonOrderingCI}
\title{Interval estimators}
\usage{
IntervalEstimator(two_sided, l1, u1, l2, u2, label)

VirtualIntervalEstimator()

RepeatedCI(two_sided = TRUE)

StagewiseCombinationFunctionOrderingCI(two_sided = TRUE)

MLEOrderingCI(two_sided = TRUE)

LikelihoodRatioOrderingCI(two_sided = TRUE)

ScoreTestOrderingCI(two_sided = TRUE)

NaiveCI(two_sided = TRUE)

NeymanPearsonOrderingCI(two_sided = TRUE, mu0 = 0, mu1 = 0.4)
}
\arguments{
\item{two_sided}{logical indicating whether the confidence interval is two-sided.}

\item{l1}{functional representation of the lower boundary of the interval in the early futility and efficacy regions.}

\item{u1}{functional representation of the upper boundary of the interval in the early futility and efficacy regions.}

\item{l2}{functional representation of the lower boundary of the interval in the continuation region.}

\item{u2}{functional representation of the upper boundary of the interval in the continuation region.}

\item{label}{name of the estimator. Used in printing methods.}

\item{mu0}{expected value of the normal distribution under the null hypothesis.}

\item{mu1}{expected value of the normal distribution under the null hypothesis.}
}
\value{
An object of class \code{IntervalEstimator}.
}
\description{
Interval estimators
}
\examples{
IntervalEstimator(
  two_sided = FALSE,
  l1 = \(smean1, ...) smean1 - 1,
  u1 = \(smean1, ...) smean1 + 1,
  l2 = \(smean2, ...) smean2 - 1,
  u2 = \(smean2, ...) smean2 + 1,
  label="My custom p-value")
}
