setGeneric("get_pdf", function(prior) standardGeneric("get_pdf"))
setGeneric("get_logpdf", function(prior) standardGeneric("get_logpdf"))
setGeneric("get_bounds", function(prior, infq) standardGeneric("get_bounds"))
setGeneric("get_mean", function(prior, infq) standardGeneric("get_mean"))
### Remove once adoptr is back on CRAN ###
setClass("Prior")
### end remove ###
setClass("NormalPrior", contains = "Prior", slots = c(mu = "numeric", sigma = "numeric"))
#' Normal prior distribution for the parameter mu
#'
#' @param mu mean of prior distribution.
#' @param sigma standard deviation of the prior distribution.
#'
#' @returns an object of class \code{NormalPrior}. This object can be supplied
#' as the argument \code{mu} of the \code{\link{evaluate_estimator}} function
#' to calculate performance scores weighted by a prior.
#'
#' @export
#'
#' @examples
#' NormalPrior(mu = 0, sigma = 1)
NormalPrior <- function(mu = 0, sigma = 1) {
 new("NormalPrior",
     mu = mu,
     sigma = sigma
     )
}
setClass("UniformPrior", contains = "Prior", slots = c(min = "numeric", max = "numeric"))
#' Uniform prior distribution for the parameter mu
#'
#' @param min minimum of support interval.
#' @param max maximum of support interval.
#'
#' @returns an object of class \code{UniformPrior}. This object can be supplied
#' as the argument \code{mu} of the \code{\link{evaluate_estimator}} function
#' to calculate performance scores weighted by a prior.
#'
#' @export
#'
#' @examples
#' UniformPrior(min = -1, max = 1)
UniformPrior <- function(min = -1, max = 1) {
  new("UniformPrior",
      min = min,
      max = max
  )
}
setMethod("get_pdf", signature = "NormalPrior",
          function(prior) function(x) dnorm(x, mean = prior@mu, sd = prior@sigma))

#' @importFrom stats dunif
setMethod("get_pdf", signature = "UniformPrior",
          function(prior) function(x) dunif(x, min = prior@min, max = prior@max))
#' @importFrom stats dnorm
setMethod("get_logpdf", signature = "NormalPrior",
          function(prior) function(x) dnorm(x, mean = prior@mu, sd = prior@sigma, log = TRUE))

setMethod("get_logpdf", signature = "UniformPrior",
          function(prior) function(x) dunif(x, min = prior@min, max = prior@max, log = TRUE))

setMethod("get_bounds", signature = "NormalPrior",
          function(prior, infq) c(qnorm(infq, mean = prior@mu, sd = prior@sigma, lower.tail = TRUE),
                                  qnorm(infq, mean = prior@mu, sd = prior@sigma, lower.tail = FALSE)))

setMethod("get_bounds", signature = "UniformPrior",
          function(prior, infq) c(prior@min, prior@max))


setMethod("get_mean", signature = "NormalPrior",
          function(prior) prior@mu)

setMethod("get_mean", signature = "UniformPrior",
          function(prior) mean(c(prior@min, prior@max)))
