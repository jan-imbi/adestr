% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adestr_package.R
\docType{package}
\name{adestr}
\alias{adestr}
\alias{adestr-package}
\title{adestr}
\description{
Point estimates, confidence intervals, and p-values for optimal adaptive two-stage designs.
}
\details{
This package implements methods to \link[adestr:evaluate_estimator]{evaluate the performance characteristics} of
various \link[adestr:PointEstimator]{point} and \link[adestr:IntervalEstimator]{interval} estimators for optimal adaptive two-stage designs.
Specifically, this package is written to interface with trial designs created by the \code{adoptr} package
\insertCite{kunzmann2021adoptr,pilz2021optimal}{adestr}.
Apart from the a priori evaluation of performance characteristics, this package also allows for the
\link[adestr:analyze]{calculation of the values of the estimators} given real datasets, and it implements methods
to calculate \link[adestr:PValue]{p-values}.
}
\references{
\insertAllCited{}
}
\seealso{
\code{\link[adestr:evaluate_estimator]{evaluate_estimator}}

\code{\link[adestr:analyze]{analyze}}

\code{\link[adestr:Statistic]{Statistic}} \code{\link[adestr:PointEstimator]{PointEstimator}} \code{\link[adestr:IntervalEstimator]{IntervalEstimator}} \code{\link[adestr:PValue]{PValue}}

\code{\link[adestr:plot,EstimatorScoreResultList-method]{plot}} \code{\link[adestr:plot_p]{plot_p}}

\url{https://jan-imbi.github.io/adestr/}
}
\author{
\strong{Maintainer}: Jan Meis \email{meis@imbi.uni-heidelberg.de} (\href{https://orcid.org/0000-0001-5407-7220}{ORCID})

Other contributors:
\itemize{
  \item Martin Maechler \email{maechler@stat.math.ethz.ch} (\href{https://orcid.org/0000-0002-8685-9910}{ORCID}) (Original author of monoSpl.c (from the 'stats' package).) [copyright holder]
}

}
