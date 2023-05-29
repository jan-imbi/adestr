setGeneric("get_pdf", function(prior) standardGeneric("get_pdf"))
setGeneric("get_logpdf", function(prior) standardGeneric("get_logpdf"))
setGeneric("get_bounds", function(prior, infq) standardGeneric("get_bounds"))
setGeneric("get_mean", function(prior, infq) standardGeneric("get_mean"))

#' @export
NormalPrior <- setClass("NormalPrior", contains = "Prior", slots = c(mu = "numeric", sigma = "numeric"))
#' @export
UniformPrior <- setClass("UniformPrior", contains = "Prior", slots = c(min = "numeric", max = "numeric"))

setMethod("get_pdf", signature = "NormalPrior",
          function(prior) function(x) dnorm(x, mean = prior@mu, sd = prior@sigma))

setMethod("get_pdf", signature = "UniformPrior",
          function(prior) function(x) dunif(x, min = prior@min, max = prior@max))

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







