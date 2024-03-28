#' Performance scores for point and interval estimators
#' @description
#' These classes encode various metrics which can be used to evaluate
#' the performance characteristics of point and interval estimators.
#'
#' @param interval confidence interval with respect to which centrality of a point
#' estimator should be evaluated.
#'
#' @details
#' # Details on the implemented estimators
#' In the following, precise definitions of the performance scores implemented
#' in \code{\link{adestr}}
#' are given. To this end,
#' let \eqn{\hat{\mu}} denote a point estimator, (\eqn{\hat{l}}, \eqn{\hat{u}})
#' an interval estimator, denote the expected value of a random variable
#' by \eqn{\mathbb{E}}, the probability of an event by \eqn{P},
#' and let \eqn{\mu} be the real value of the underlying
#' parameter to be estimated.
#'
#' ## Scores for point estimators (\code{PointEstimatorScore}):
#' * \code{Expectation()}: \eqn{\mathbb{E}[\hat{\mu}]}
#' * \code{Bias()}: \eqn{\mathbb{E}[\hat{\mu} - \mu]}
#' * \code{Variance()}: \eqn{\mathbb{E}[(\hat{\mu} - \mathbb{E}[\hat{\mu}])^2]}
#' * \code{MSE()}: \eqn{\mathbb{E}[(\hat{\mu} - mu)^2]}
#' * \code{OverestimationProbability()}: \eqn{P(\hat{\mu} > \mu)}
#' * \code{Centrality(interval)}: \eqn{\mathbb{E}[(\hat{\mu} - \hat{l}) + (\hat{\mu} - \hat{u}]}
#' ## Scores for confidence intervals (\code{IntervalEstimatorScore}):
#' * \code{Coverage()}: \eqn{P(\hat{l} \leq \mu \leq \hat{u})}
#' * \code{Width()}: \eqn{\mathbb{E}[\hat{u} - \hat{l}]}
#' * \code{TestAgreement()}: \eqn{P\left( \left(\{0 < \hat{l} \text{ and } (c_{1, e} < Z_1 \text{ or } c_{2}(Z_1) < Z_2 ) \right) \text{ or } \left(\{\hat{l} \leq 0  \text{ and } ( Z_1 < c_{1, f} \text{ or } Z_2 \leq c_{2}(Z_1))\}\right)\right)}
#'
#' @md
#' @slot label name of the performance score. Used in printing methods.
#'
#' @returns an object of class \code{EstimatorScore}. This class signals that
#' an object can be used with the \code{\link{evaluate_estimator}} function.
#' @export
#' @aliases EstimatorScore
#' @seealso \code{\link{evaluate_estimator}}
#'
#' @inherit evaluate_estimator examples
setClass("EstimatorScore", slots = c(label = "character"))
setClass("PointEstimatorScore", contains = "EstimatorScore")
setClass("IntervalEstimatorScore", contains = "EstimatorScore")
EstimatorScoreResult <- setClass("EstimatorScoreResult", slots = c(score = "list", estimator = "Statistic", data_distribution = "DataDistribution",
                                                                   design = "TwoStageDesign", mu = "ANY", sigma = "numeric",
                                                                   results = "list", integrals = "list"))
setClass("EstimatorScoreResultList", contains = "list")
EstimatorScoreResultList <- function(...) {
  r <- list(...)
  class(r) <- c("EstimatorScoreResultList", "list")
  r
}

#' Combine EstimatoreScoreResult objects into a list
#'
#' Creates an object of class EstimatoreScoreResultList,
#' which is a basically list with the respective
#' EstimatoreScoreResult objects.
#'
#' @param x an object of class EstimatorScoreResult.
#' @param ... additional arguments passed along to the \code{\link{list}} function
#' @return an object of class EstimatoreScoreResultList.
setMethod("c", signature("EstimatorScoreResult"), definition =
            function(x, ...) {
              EstimatorScoreResultList(x, ...)
            })
#' @inherit c,EstimatorScoreResult-method
setMethod("c", signature("EstimatorScoreResultList"), definition =
            function(x, ...) {
              EstimatorScoreResultList(x, ...)
            })

#' Evaluate performance characteristics of an estimator
#'
#' This function evaluates an \code{\link{EstimatorScore}} for a \code{\link{PointEstimator}}
#' or and \code{\link{IntervalEstimator}} by integrating over the sampling distribution.
#'
#' @details
#' ## General
#'
#' First, a functional representation of the integrand is created by combining information
#' from the \code{\link{EstimatorScore}} object (\code{score}) and the \code{\link{PointEstimator}} or
#' \code{\link{IntervalEstimator}} object (\code{estimator}).
#' The sampling distribution of a design is determined by the \code{TwoStageDesign} object
#' (\code{design}) and the \code{DataDistribution} object (\code{data_distribution}),
#' as well as the assumed parameters \eqn{\mu} (mu) and \eqn{\sigma} (sigma).
#' The other parameters control various details of the integration problem.
#'
#' ## Other parameters
#'
#' For a two-armed \code{data_distribution},
#' if \code{use_full_twoarm_sampling_distribution} is \code{TRUE}, the sample means
#' for both groups are integrated independently. If \code{use_full_twoarm_sampling_distribution}
#' is \code{FALSE}, only the difference in sample means is integrated.
#'
#' \code{true_parameter} controls which parameters is supposed to be estimated. This
#' is usually \code{mu}, but could be set to \code{sigma} if one is interested in
#' estimating the standard deviation.
#'
#' If the parameter \code{exact} is set to \code{FALSE}
#' (the default), the continuous version of the second-stage sample-size function \code{n2}
#' is used. Otherwise, an integer valued version of that function will be used,
#' though this is considerably slower.
#'
#' The parameters \code{early_futility_part},
#' \code{continuation_part} and \code{early_efficacy_part} control which parts
#' of the sample-space should be integrated over (all default to \code{TRUE}).
#' They can be used in conjunction with the parameter \code{conditional_integral},
#' which enables the calculation of the expected value of performance score conditional
#' on reaching any of the selected integration regions.
#'
#' Lastly, the paramters
#' \code{tol}, \code{maxEval}, and \code{absError} control the integration accuracy.
#' They are handed down to the \code{\link{hcubature}} function.
#' @md
#'
#' @param score performance measure to evaluate.
#' @param estimator object of class \code{PointEstimator}, \code{IntervalEstimator} or \code{PValue}.
#' @param data_distribution object of class \code{Normal} or \code{Student}.
#' @param use_full_twoarm_sampling_distribution logical indicating whether this estimator is intended to be used
#' with the full sampling distribution in a two-armed trial.
#' @param design object of class \code{TwoStageDesign}.
#' @param true_parameter true value of the parameter (used e.g. when evaluating bias).
#' @param mu expected value of the underlying normal distribution.
#' @param sigma assumed standard deviation.
#' @param tol relative tolerance.
#' @param maxEval maximum number of iterations.
#' @param absError absolute tolerance.
#' @param exact logical indicating usage of exact n2 function.
#' @param early_futility_part include early futility part of integral.
#' @param continuation_part include continuation part of integral.
#' @param early_efficacy_part include early efficacy part of integral.
#' @param conditional_integral treat integral as a conditional integral.
#'
#' @return an object of class \code{EstimatorScoreResult}
#' containing the values of the evaluated \code{\link{EstimatorScore}} and
#' information about the setting for which they were calculated
#' (e.g. the \code{estimator}, \code{data_distribution}, \code{design}, \code{mu}, and \code{sigma}).
#'
#' @seealso [EstimatorScore]
#' @seealso [PointEstimator] [IntervalEstimator]
#' @seealso \link[adestr:plot,EstimatorScoreResultList-method]{plot}
#' @export
#'
#' @examples
#' evaluate_estimator(
#'   score = MSE(),
#'   estimator = SampleMean(),
#'   data_distribution = Normal(FALSE),
#'   design = get_example_design(),
#'   mu = c(0, 0.3, 0.6),
#'   sigma = 1,
#'   exact = FALSE
#' )
#'
#' evaluate_estimator(
#'   score = Coverage(),
#'   estimator = StagewiseCombinationFunctionOrderingCI(),
#'   data_distribution = Normal(FALSE),
#'   design = get_example_design(),
#'   mu = c(0, 0.3),
#'   sigma = 1,
#'   exact = FALSE
#' )
#'
setGeneric("evaluate_estimator", function(score,
                                          estimator,
                                          data_distribution,
                                          use_full_twoarm_sampling_distribution = FALSE,
                                          design,
                                          true_parameter = mu,
                                          mu,
                                          sigma,
                                          tol = getOption("adestr_tol_outer", default = .adestr_options[["adestr_tol_outer"]]),
                                          maxEval = getOption("adestr_maxEval_outer", default = .adestr_options[["adestr_maxEval_outer"]]),
                                          absError = getOption("adestr_absError_outer", default = .adestr_options[["adestr_absError_outer"]]),
                                          exact = FALSE,
                                          early_futility_part = TRUE,
                                          continuation_part = TRUE,
                                          early_efficacy_part = TRUE,
                                          conditional_integral = FALSE) standardGeneric("evaluate_estimator"))

#' @inherit evaluate_estimator
#' @name evaluate_estimator-methods
NULL

#' @rdname evaluate_estimator-methods
setMethod("evaluate_estimator", signature("PointEstimatorScore", "IntervalEstimator"),
          function(score,
                   estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   true_parameter,
                   mu,
                   sigma,
                   tol,
                   maxEval,
                   absError,
                   exact,
                   early_futility_part,
                   continuation_part,
                   early_efficacy_part,
                   conditional_integral) {
            stop("Cannot evaluate PointEstimatorScore on IntervalEstimator.")
          })
#' @rdname evaluate_estimator-methods
setMethod("evaluate_estimator", signature("IntervalEstimatorScore", "PointEstimator"),
          function(score,
                   estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   true_parameter,
                   mu,
                   sigma,
                   tol,
                   maxEval,
                   absError,
                   exact,
                   early_futility_part,
                   continuation_part,
                   early_efficacy_part,
                   conditional_integral) {
            stop("Cannot evaluate IntervalEstimatorScore on PointEstimator.")
          })
#' @rdname evaluate_estimator-methods
setMethod("evaluate_estimator", signature("list", "Estimator"),
          function(score,
                   estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   true_parameter,
                   mu,
                   sigma,
                   tol,
                   maxEval,
                   absError,
                   exact,
                   early_futility_part,
                   continuation_part,
                   early_efficacy_part,
                   conditional_integral) {
            .results <- lapply(score, evaluate_estimator,
                               estimator = estimator,
                               data_distribution = data_distribution,
                               use_full_twoarm_sampling_distribution = use_full_twoarm_sampling_distribution,
                               design = design,
                               true_parameter = true_parameter,
                               mu = mu,
                               sigma = sigma,
                               tol = tol,
                               maxEval = maxEval,
                               absError = absError,
                               exact = exact,
                               early_futility_part = early_futility_part,
                               continuation_part = continuation_part,
                               early_efficacy_part = early_efficacy_part,
                               conditional_integral = conditional_integral)
            results <- do.call("c", lapply(.results, \(x) x@results))
            integrals <-  do.call("c", lapply(.results, \(x) x@integrals))
            return(
              EstimatorScoreResult(
                score = score,
                estimator = estimator,
                data_distribution = data_distribution,
                design = design,
                mu = mu,
                sigma = sigma,
                results = results,
                integrals = integrals
              )
            )
          })

.evaluate_estimator <- function(score,
                                estimator,
                                data_distribution,
                                use_full_twoarm_sampling_distribution,
                                design,
                                generate_g1,
                                generate_g2,
                                true_parameter,
                                mu,
                                sigma,
                                tol,
                                maxEval,
                                absError,
                                exact,
                                early_futility_part,
                                continuation_part,
                                early_efficacy_part,
                                conditional_integral){
  two_armed <- data_distribution@two_armed
  integrals <- lapply(
    seq_along(mu),
    \(i) {
      m <- mu[[i]]
      if (length(true_parameter) > 1)
        tp <- true_parameter[[i]]
      else
        tp <- true_parameter
      g1 <- generate_g1(tp)
      g2 <- generate_g2(tp)
      integrate_over_sample_space(
        data_distribution,
        use_full_twoarm_sampling_distribution,
        design,
        g1,
        g2,
        m,
        sigma,
        tol,
        maxEval,
        absError,
        exact,
        early_futility_part,
        continuation_part,
        early_efficacy_part,
        conditional_integral
      )
    })
  results <- list(sapply(integrals, \(x)x$overall_integral$integral))
  names(results) <- score@label
  integrals <- integrals
  names(integrals) <- score@label
  return(
    EstimatorScoreResult(
      score = list(score),
      estimator = estimator,
      data_distribution = data_distribution,
      design = design,
      mu = mu,
      sigma = sigma,
      results = results,
      integrals = integrals
    )
  )
}

.evaluate_estimator_prior <- function(score,
                                estimator,
                                data_distribution,
                                use_full_twoarm_sampling_distribution,
                                design,
                                g1,
                                g2,
                                true_parameter,
                                mu,
                                sigma,
                                tol,
                                maxEval,
                                absError,
                                exact,
                                early_futility_part,
                                continuation_part,
                                early_efficacy_part,
                                conditional_integral){
  two_armed <- data_distribution@two_armed
  integrals <- list(
    integrate_over_sample_space(
      data_distribution,
      use_full_twoarm_sampling_distribution,
      design,
      g1,
      g2,
      mu,
      sigma,
      tol,
      maxEval,
      absError,
      exact,
      early_futility_part,
      continuation_part,
      early_efficacy_part,
      conditional_integral
    )
  )
  results <- list(sapply(integrals, \(x)x$overall_integral$integral))
  names(results) <- score@label
  integrals <- integrals
  names(integrals) <- score@label
  return(
    EstimatorScoreResult(
      score = list(score),
      estimator = estimator,
      data_distribution = data_distribution,
      design = design,
      mu = mu,
      sigma = sigma,
      results = results,
      integrals = integrals
    )
  )
}

setClass("Expectation", contains = "PointEstimatorScore")
#' @rdname EstimatorScore-class
#' @export
Expectation <- function() new("Expectation", label = "Expectation")
#' @rdname evaluate_estimator-methods
setMethod("evaluate_estimator", signature("Expectation", "PointEstimator"),
          function(score,
                   estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   true_parameter,
                   mu,
                   sigma,
                   tol,
                   maxEval,
                   absError,
                   exact,
                   early_futility_part,
                   continuation_part,
                   early_efficacy_part,
                   conditional_integral) {
            stagewise_estimators <- get_stagewise_estimators(estimator = estimator,
                                                             data_distribution =  data_distribution,
                                                             use_full_twoarm_sampling_distribution = use_full_twoarm_sampling_distribution,
                                                             design = design, sigma = sigma, exact = exact)
            g1 <- stagewise_estimators[[1L]]
            g2 <- stagewise_estimators[[2L]]
            .evaluate_estimator(
              score,
              estimator,
              data_distribution,
              use_full_twoarm_sampling_distribution,
              design,
              generate_g1 = \(tp) g1,
              generate_g2 = \(tp) g2,
              true_parameter,
              mu,
              sigma,
              tol,
              maxEval,
              absError,
              exact,
              early_futility_part,
              continuation_part,
              early_efficacy_part,
              conditional_integral)
          })


setClass("Bias", contains = "PointEstimatorScore")
#' @rdname EstimatorScore-class
#' @export
Bias <- function() new("Bias", label = "Bias")
#' @rdname evaluate_estimator-methods
setMethod("evaluate_estimator", signature("Bias", "PointEstimator"),
          function(score,
                   estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   true_parameter,
                   mu,
                   sigma,
                   tol,
                   maxEval,
                   absError,
                   exact,
                   early_futility_part,
                   continuation_part,
                   early_efficacy_part,
                   conditional_integral) {
            res <- evaluate_estimator(Expectation(),
                                      estimator,
                                      data_distribution,
                                      use_full_twoarm_sampling_distribution,
                                      design,
                                      true_parameter,
                                      mu,
                                      sigma,
                                      tol,
                                      maxEval,
                                      absError,
                                      exact,
                                      early_futility_part,
                                      continuation_part,
                                      early_efficacy_part,
                                      conditional_integral)

            res@results$Bias <- res@results$Expectation - true_parameter
            res@integrals$Bias <- res@integrals$Expectation
            res@score <- list(Bias())
            return(res)
          })

setClass("Variance", contains = "PointEstimatorScore")
#' @rdname EstimatorScore-class
#' @export
Variance <- function() new("Variance", label = "Variance")
#' @rdname evaluate_estimator-methods
setMethod("evaluate_estimator", signature("Variance", "PointEstimator"),
          function(score,
                   estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   true_parameter,
                   mu,
                   sigma,
                   tol,
                   maxEval,
                   absError,
                   exact,
                   early_futility_part,
                   continuation_part,
                   early_efficacy_part,
                   conditional_integral) {
            res <- evaluate_estimator(Expectation(),
                                      estimator,
                                      data_distribution,
                                      use_full_twoarm_sampling_distribution,
                                      design,
                                      true_parameter,
                                      mu,
                                      sigma,
                                      tol,
                                      maxEval,
                                      absError,
                                      exact,
                                      early_futility_part,
                                      continuation_part,
                                      early_efficacy_part,
                                      conditional_integral)
            stagewise_estimators <- get_stagewise_estimators(estimator = estimator,
                                                             data_distribution =  data_distribution,
                                                             use_full_twoarm_sampling_distribution = use_full_twoarm_sampling_distribution,
                                                             design = design, sigma = sigma, exact = exact)
            g1 <- stagewise_estimators[[1L]]
            g2 <- stagewise_estimators[[2L]]
            res2 <- .evaluate_estimator(
              score,
              estimator,
              data_distribution,
              use_full_twoarm_sampling_distribution,
              design,
              generate_g1 = \(tp) \(...)(g1(...)-tp)^2,
              generate_g2 = \(tp) \(...)(g2(...)-tp)^2,
              true_parameter,
              mu,
              sigma,
              tol,
              maxEval,
              absError,
              exact,
              early_futility_part,
              continuation_part,
              early_efficacy_part,
              conditional_integral)
            res@results$Variance <- res2@results$Variance
            res@integrals$Variance <- res2@integrals$Variance
            res@score <- list(Variance())
            return(res)
          })
setClass("MSE", contains = "PointEstimatorScore")
#' @rdname EstimatorScore-class
#' @export
MSE <- function() new("MSE", label = "MSE")
#' @rdname evaluate_estimator-methods
setMethod("evaluate_estimator", signature("MSE", "PointEstimator"),
          function(score,
                   estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   true_parameter,
                   mu,
                   sigma,
                   tol,
                   maxEval,
                   absError,
                   exact,
                   early_futility_part,
                   continuation_part,
                   early_efficacy_part,
                   conditional_integral) {
            if (is(mu, "Prior")) {
              stagewise_estimators <- get_stagewise_estimators(estimator = estimator,
                                                               data_distribution =  data_distribution,
                                                               use_full_twoarm_sampling_distribution = use_full_twoarm_sampling_distribution,
                                                               design = design, sigma = sigma, exact = exact)
              g1 <- stagewise_estimators[[1L]]
              g2 <- stagewise_estimators[[2L]]
              res <- .evaluate_estimator_prior(
                score,
                estimator,
                data_distribution,
                use_full_twoarm_sampling_distribution,
                design,
                g1 = \(mu, ...)(g1(...)-mu)^2,
                g2 = \(mu, ...)(g2(...)-mu)^2,
                true_parameter,
                mu,
                sigma,
                tol,
                maxEval,
                absError,
                exact,
                early_futility_part,
                continuation_part,
                early_efficacy_part,
                conditional_integral)
              return(res)
            } else {
              res <- evaluate_estimator(Variance(),
                                        estimator,
                                        data_distribution,
                                        use_full_twoarm_sampling_distribution,
                                        design,
                                        true_parameter,
                                        mu,
                                        sigma,
                                        tol,
                                        maxEval,
                                        absError,
                                        exact,
                                        early_futility_part,
                                        continuation_part,
                                        early_efficacy_part,
                                        conditional_integral)
              res@results$Bias <- res@results$Expectation - true_parameter
              res@results$MSE <- res@results$Bias^2 + res@results$Variance
              res@results <- res@results[c(1,3,2,4)]
              res@integrals$Bias <- res@integrals$Expectation
              res@integrals$MSE <- lapply(names(res@integrals$Expectation), \(x){
                lapply(names(res@integrals$Expectation[[x]]), \(y) res@integrals$Expectation[[x]][[y]] + res@integrals$Variance[[x]][[y]])
              })
              names(res@integrals$MSE) <- names(res@integrals$Expectation)
              for (i in seq_along(res@integrals$MSE))
                names(res@integrals$MSE[[i]]) <- names(res@integrals$Variance[[i]])
              res@integrals <- res@integrals[c(1,3,2,4)]
              res@score <- list(MSE())
              return(res)
            }
          })
setClass("OverestimationProbability", contains = "PointEstimatorScore")
#' @rdname EstimatorScore-class
#' @export
OverestimationProbability <- function() new("OverestimationProbability", label = "Probability of overestimation")
#' @rdname evaluate_estimator-methods
setMethod("evaluate_estimator", signature("OverestimationProbability", "PointEstimator"),
          function(score,
                   estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   true_parameter,
                   mu,
                   sigma,
                   tol,
                   maxEval,
                   absError,
                   exact,
                   early_futility_part,
                   continuation_part,
                   early_efficacy_part,
                   conditional_integral) {
            stagewise_estimators <- get_stagewise_estimators(estimator = estimator,
                                                             data_distribution =  data_distribution,
                                                             use_full_twoarm_sampling_distribution = use_full_twoarm_sampling_distribution,
                                                             design = design, sigma = sigma, exact = exact)
            g1 <- stagewise_estimators[[1L]]
            g2 <- stagewise_estimators[[2L]]
            .evaluate_estimator(
              score,
              estimator,
              data_distribution,
              use_full_twoarm_sampling_distribution,
              design,
              generate_g1 = \(tp) \(...) (g1(...) > tp),
              generate_g2 = \(tp) \(...) (g2(...) > tp),
              true_parameter,
              mu,
              sigma,
              tol,
              maxEval,
              absError,
              exact,
              early_futility_part,
              continuation_part,
              early_efficacy_part,
              conditional_integral)
          })
setClass("Coverage", contains = "IntervalEstimatorScore")
#' @rdname EstimatorScore-class
#' @export
Coverage <- function() new("Coverage", label = "Coverage")
#' @rdname evaluate_estimator-methods
setMethod("evaluate_estimator", signature("Coverage", "IntervalEstimator"),
          function(score,
                   estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   true_parameter,
                   mu,
                   sigma,
                   tol,
                   maxEval,
                   absError,
                   exact,
                   early_futility_part,
                   continuation_part,
                   early_efficacy_part,
                   conditional_integral) {
            stagewise_estimators <- get_stagewise_estimators(estimator = estimator,
                                                             data_distribution =  data_distribution,
                                                             use_full_twoarm_sampling_distribution = use_full_twoarm_sampling_distribution,
                                                             design = design, sigma = sigma, exact = exact)
            l1 <- stagewise_estimators[[1L]]
            u1 <- stagewise_estimators[[2L]]
            l2 <- stagewise_estimators[[3L]]
            u2 <- stagewise_estimators[[4L]]
            .evaluate_estimator(
              score,
              estimator,
              data_distribution,
              use_full_twoarm_sampling_distribution,
              design,
              generate_g1 = \(tp) \(...) l1(...) <= tp & tp <= u1(...),
              generate_g2 = \(tp) \(...) l2(...) <= tp & tp <= u2(...),
              true_parameter,
              mu,
              sigma,
              tol,
              maxEval,
              absError,
              exact,
              early_futility_part,
              continuation_part,
              early_efficacy_part,
              conditional_integral)
          })
setClass("SoftCoverage", contains = "IntervalEstimatorScore", slots = c(shrinkage = "numeric"))
#' @rdname EstimatorScore-class
#' @param shrinkage shrinkage factor for bump function.
#' @export
SoftCoverage <- function(shrinkage = 1) new("SoftCoverage", shrinkage = shrinkage, label = "Coverage")
#' @rdname evaluate_estimator-methods
setMethod("evaluate_estimator", signature("SoftCoverage", "IntervalEstimator"),
          function(score,
                   estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   true_parameter,
                   mu,
                   sigma,
                   tol,
                   maxEval,
                   absError,
                   exact,
                   early_futility_part,
                   continuation_part,
                   early_efficacy_part,
                   conditional_integral) {
            stagewise_estimators <- get_stagewise_estimators(estimator = estimator,
                                                             data_distribution =  data_distribution,
                                                             use_full_twoarm_sampling_distribution = use_full_twoarm_sampling_distribution,
                                                             design = design, sigma = sigma, exact = exact)
            l1 <- stagewise_estimators[[1L]]
            u1 <- stagewise_estimators[[2L]]
            l2 <- stagewise_estimators[[3L]]
            u2 <- stagewise_estimators[[4L]]
            .evaluate_estimator(
              score,
              estimator,
              data_distribution,
              use_full_twoarm_sampling_distribution,
              design,
              generate_g1 = \(tp) \(...) {
                ll <- l1(...)
                uu <- u1(...)
                pmax(ll <= tp & tp <= uu, exp(-score@shrinkage * (ll - tp)^2), exp(-score@shrinkage * (uu - tp)^2))
              },
              generate_g2 = \(tp) \(...) l2(...) <= tp & tp <= u2(...),
              true_parameter,
              mu,
              sigma,
              tol,
              maxEval,
              absError,
              exact,
              early_futility_part,
              continuation_part,
              early_efficacy_part,
              conditional_integral)
          })



setClass("Width", contains = "IntervalEstimatorScore")
#' @rdname EstimatorScore-class
#' @export
Width <- function() new("Width", label = "Width")
#' @rdname evaluate_estimator-methods
setMethod("evaluate_estimator", signature("Width", "IntervalEstimator"),
          function(score,
                   estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   true_parameter,
                   mu,
                   sigma,
                   tol,
                   maxEval,
                   absError,
                   exact,
                   early_futility_part,
                   continuation_part,
                   early_efficacy_part,
                   conditional_integral) {
            stagewise_estimators <- get_stagewise_estimators(estimator = estimator,
                                                             data_distribution =  data_distribution,
                                                             use_full_twoarm_sampling_distribution = use_full_twoarm_sampling_distribution,
                                                             design = design, sigma = sigma, exact = exact)
            l1 <- stagewise_estimators[[1L]]
            u1 <- stagewise_estimators[[2L]]
            l2 <- stagewise_estimators[[3L]]
            u2 <- stagewise_estimators[[4L]]
            .evaluate_estimator(
              score,
              estimator,
              data_distribution,
              use_full_twoarm_sampling_distribution,
              design,
              generate_g1 = \(tp) \(...) abs(u1(...) - l1(...)),
              generate_g2 = \(tp) \(...) abs(u2(...) - l2(...)),
              true_parameter,
              mu,
              sigma,
              tol,
              maxEval,
              absError,
              exact,
              early_futility_part,
              continuation_part,
              early_efficacy_part,
              conditional_integral)
          })
setClass("TestAgreement", contains = "IntervalEstimatorScore")
#' @rdname EstimatorScore-class
#' @export
TestAgreement <- function() new("TestAgreement", label = "Agreement with test decision")
#' @rdname evaluate_estimator-methods
setMethod("evaluate_estimator", signature("TestAgreement", "IntervalEstimator"),
          function(score,
                   estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   true_parameter,
                   mu,
                   sigma,
                   tol,
                   maxEval,
                   absError,
                   exact,
                   early_futility_part,
                   continuation_part,
                   early_efficacy_part,
                   conditional_integral) {
            design <- TwoStageDesignWithCache(design)
            stagewise_estimators <- get_stagewise_estimators(estimator = estimator,
                                                             data_distribution =  data_distribution,
                                                             use_full_twoarm_sampling_distribution = use_full_twoarm_sampling_distribution,
                                                             design = design, sigma = sigma, exact = exact)
            l1 <- stagewise_estimators[[1L]]
            u1 <- stagewise_estimators[[2L]]
            l2 <- stagewise_estimators[[3L]]
            u2 <- stagewise_estimators[[4L]]

            generate_g1 <- \(tp) \(design, smean1, n1, sigma, two_armed, ...) {
              stop("Unreachable.")
            }
            generate_g2 <- \(tp) \(design, smean1, smean2, n1, n2, sigma, two_armed, ...) {
              stop("Unreachable.")
            }
            if (is(data_distribution, "Student")) {
              generate_g1 = \(tp) \(design, smean1, svar1, n1, two_armed, ...) {
                t1 <- smean_to_t(smean1, svar1, n1, two_armed)
                !xor(design@c1e < t1, 0 < l1(design = design, smean1 = smean1, svar1 = svar1, n1 = n1, two_armed = two_armed, ...))
              }
              generate_g2 <- \(tp) \(design, smean1, svar1, smean2, svar2, n1, n2, two_armed, ...) {
                t1 <- smean_to_t(smean1, svar1, n1, two_armed)
                t2 <- smean_to_t(smean2, svar2, n2, two_armed)
                c2 <- c2_extrapol(design, t1)
                !xor(c2 < t2, 0 < l2(design = design, smean1 = smean1, svar1 = svar1, smean2 = smean2, svar2 = svar2, n1 = n1, n2 = n2, two_armed = two_armed, ...))
              }
            } else if (is(data_distribution, "Normal")) {
              generate_g1 = \(tp) \(design, smean1, n1, sigma, two_armed, ...) {
                z1 <- smean_to_z(smean1, n1, sigma, two_armed)
                !xor(design@c1e < z1, 0 < l1(design = design, smean1 = smean1, n1 = n1, sigma = sigma, two_armed = two_armed, ...))
              }
              generate_g2 <- \(tp) \(design, smean1, smean2, n1, n2, sigma, two_armed, ...) {
                z1 <- smean_to_z(smean1, n1, sigma, two_armed)
                z2 <- smean_to_z(smean2, n2, sigma, two_armed)
                c2 <- c2_extrapol(design, z1)
                !xor(c2 < z2, 0 < l2(design = design, smean1 = smean1, smean2 = smean2, n1 = n1, n2 = n2, sigma = sigma, two_armed = two_armed, ...))
              }
            } else {
              stop("This data_distribution class is not supported.")
            }
            .evaluate_estimator(
              score,
              estimator,
              data_distribution,
              use_full_twoarm_sampling_distribution,
              design,
              generate_g1,
              generate_g2,
              true_parameter,
              mu,
              sigma,
              tol,
              maxEval,
              absError,
              exact,
              early_futility_part,
              continuation_part,
              early_efficacy_part,
              conditional_integral)
          })

#' @rdname evaluate_estimator-methods
setMethod("evaluate_estimator", signature("TestAgreement", "PValue"),
          function(score,
                   estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   true_parameter,
                   mu,
                   sigma,
                   tol,
                   maxEval,
                   absError,
                   exact,
                   early_futility_part,
                   continuation_part,
                   early_efficacy_part,
                   conditional_integral) {
            if (missing(true_parameter))
              true_parameter <- 0.025
            design <- TwoStageDesignWithCache(design)
            stagewise_estimators <- get_stagewise_estimators(estimator = estimator,
                                                             data_distribution =  data_distribution,
                                                             use_full_twoarm_sampling_distribution = use_full_twoarm_sampling_distribution,
                                                             design = design, sigma = sigma, exact = exact)
            g1 <- stagewise_estimators[[1L]]
            g2 <- stagewise_estimators[[2L]]

            if (is(data_distribution, "Student")) {
              generate_g1 = \(tp) \(design, smean1, svar1, n1, two_armed, ...) {
                t1 <- smean_to_t(smean1, svar1, n1, two_armed)
                !xor(design@c1e < t1, g1(design = design, smean1 = smean1, svar1 = svar1, n1 = n1, two_armed = two_armed, ...) < true_parameter)
              }
              generate_g2 <- \(tp) \(design, smean1, svar1, smean2, svar2, n1, n2, two_armed, ...) {
                t1 <- smean_to_t(smean1, svar1, n1, two_armed)
                t2 <- smean_to_t(smean2, svar2, n2, two_armed)
                c2 <- c2_extrapol(design, t1)
                !xor(c2 < t2, g2(design = design, smean1 = smean1, svar1 = svar1, smean2 = smean2, svar2 = svar2, n1 = n1, n2 = n2, two_armed = two_armed, ...) < true_parameter)
              }
            } else if (is(data_distribution, "Normal")) {
              generate_g1 = \(tp) \(design, smean1, n1, sigma, two_armed, ...) {
                z1 <- smean_to_z(smean1, n1, sigma, two_armed)
                !xor(design@c1e < z1, g1(design = design, smean1 = smean1, n1 = n1, sigma = sigma, two_armed = two_armed, ...) < true_parameter)
              }
              generate_g2 <- \(tp) \(design, smean1, smean2, n1, n2, sigma, two_armed, ...) {
                z1 <- smean_to_z(smean1, n1, sigma, two_armed)
                z2 <- smean_to_z(smean2, n2, sigma, two_armed)
                c2 <- c2_extrapol(design, z1)
                !xor(c2 < z2, g2(design = design, smean1 = smean1, smean2 = smean2, n1 = n1, n2 = n2, sigma = sigma, two_armed = two_armed, ...) < true_parameter)
              }
            } else {
              stop("This data_distribution class is not supported.")
            }
            .evaluate_estimator(
              score,
              estimator,
              data_distribution,
              use_full_twoarm_sampling_distribution,
              design,
              generate_g1,
              generate_g2,
              true_parameter,
              mu,
              sigma,
              tol,
              maxEval,
              absError,
              exact,
              early_futility_part,
              continuation_part,
              early_efficacy_part,
              conditional_integral)
          })


setClass("Centrality", contains = "PointEstimatorScore",
         slots = list(interval = "IntervalEstimator"))
#' @rdname EstimatorScore-class
#' @export
Centrality <- function(interval = NULL) new("Centrality", label = sprintf("Centrality with respect to %s", toString(interval)), interval = interval)
#' @rdname evaluate_estimator-methods
setMethod("evaluate_estimator", signature("Centrality", "PointEstimator"),
          function(score,
                   estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   true_parameter,
                   mu,
                   sigma,
                   tol,
                   maxEval,
                   absError,
                   exact,
                   early_futility_part,
                   continuation_part,
                   early_efficacy_part,
                   conditional_integral) {
            design <- TwoStageDesignWithCache(design)
            stagewise_estimators <- get_stagewise_estimators(estimator = estimator,
                                                             data_distribution =  data_distribution,
                                                             use_full_twoarm_sampling_distribution = use_full_twoarm_sampling_distribution,
                                                             design = design, sigma = sigma, exact = exact)
            stagewise_intervals <- get_stagewise_estimators(estimator = score@interval,
                                                             data_distribution =  data_distribution,
                                                             use_full_twoarm_sampling_distribution = use_full_twoarm_sampling_distribution,
                                                             design = design, sigma = sigma, exact = exact)
            g1 <- stagewise_estimators[[1L]]
            g2 <- stagewise_estimators[[2L]]
            l1 <- stagewise_intervals[[1L]]
            u1 <- stagewise_intervals[[2L]]
            l2 <- stagewise_intervals[[3L]]
            u2 <- stagewise_intervals[[4L]]
            generate_g1 = \(tp) \(...) ((g1(...) - l1(...)) + (g1(...) - u1(...)))
            generate_g2 = \(tp) \(...) ((g2(...) - l2(...)) + (g2(...) - u2(...)))
            .evaluate_estimator(
              score,
              estimator,
              data_distribution,
              use_full_twoarm_sampling_distribution,
              design,
              generate_g1,
              generate_g2,
              true_parameter,
              mu,
              sigma,
              tol,
              maxEval,
              absError,
              exact,
              early_futility_part,
              continuation_part,
              early_efficacy_part,
              conditional_integral)
          })

#' Evaluate different scenarios in parallel
#'
#' This function takes a list of lists of scores, a list of lists of estimators,
#' and lists lists of various other design parameters. Each possible combination
#' of the elements of the respective sublists is then used to create separate
#' scenarios.
#' These scenarios are than evaluated independelty of each other,
#' allowing for parallelization via the \code{\link[future:future]{future}} framework. For each scenario,
#' one call to the \code{\link{evaluate_estimator}} function is made.
#'
#' Concretely, the cross product of
#' the first sublist of scores and the first sublist of estimators and the other parameters
#' is calculated. Then the cross product of the second sublist of scores, estimators
#' and other design parameters is calculated. All of these cross products together
#' make up the set of all scenarios. The combinations say the first sublist of scores
#' and the second sublist of estimators are not considered.
#'
#' @param score_lists a list of lists of estimator scores.
#' @param estimator_lists a list of lists of estimators.
#' @param data_distribution_lists a list of lists of data distributions.
#' @param use_full_twoarm_sampling_distribution_lists a list of lists of use_full_twoarm_sampling_distribution_lists parameters.
#' @param design_lists a list of lists of designs.
#' @param true_parameter_lists a list of lists of true parameters.
#' @param mu_lists a list of lists of mu vectors.
#' @param sigma_lists a list of lists of sigma values.
#' @param tol_lists a list of lists of relative tolerances.
#' @param maxEval_lists a list of lists of maxEval boundaries.
#' @param absError_lists a list of lists of absError boundaries.
#' @param exact_lists a list of lists of `exact` parameters.
#' @param early_futility_part_lists a list of lists of `early_futility_part_lists` parameters.
#' @param continuation_part_lists a list of lists of `continuation_part_lists` parameters.
#' @param early_efficacy_part_lists a list of lists of `early_efficacy_part_lists` parameters.
#' @param conditional_integral_lists a list of lists of `conditional_integral_lists` parameters.
#'
#' @returns a list of data.frames containing the results for the respective scenarios.
#' @export
#'
#' @examples
#' res <-evaluate_scenarios_parallel(
#'  score_lists = list(c(MSE(), OverestimationProbability())),
#'  estimator_lists =  list(c(SampleMean(), FirstStageSampleMean())),
#'  data_distribution_lists = list(c(Normal(FALSE), Normal(TRUE))),
#'  design_lists =  list(c(get_example_design())),
#'  mu_lists = list(c(-1, 0, 1)),
#'  sigma_lists = list(1)
#' )
#'
#' @importFrom future.apply future_apply
#' @importFrom progressr progressor
#' @seealso [evaluate_estimator]
evaluate_scenarios_parallel <- function(score_lists,
                                        estimator_lists,
                                        data_distribution_lists,
                                        use_full_twoarm_sampling_distribution_lists,
                                        design_lists,
                                        true_parameter_lists,
                                        mu_lists,
                                        sigma_lists,
                                        tol_lists,
                                        maxEval_lists,
                                        absError_lists,
                                        exact_lists,
                                        early_futility_part_lists,
                                        continuation_part_lists,
                                        early_efficacy_part_lists,
                                        conditional_integral_lists) {
  passed <- names(as.list(match.call())[-1])
  argnames_wo_list <- sapply(passed, \(x)substr(x, 1, nchar(x)-nchar("_lists")))
  input_args <- mget(passed)
  for (i in seq_along(input_args)) {
    if (!is.list(input_args[[i]]) || (!is.list(input_args[[i]][[1L]]) && !is.vector(input_args[[i]][[1L]])))
      input_args[[i]] <- list(input_args[[i]])
    for (j in seq_along(input_args[[i]])) {
      if (!is.list(input_args[[i]][[j]]) && !is.vector(input_args[[i]][[j]])) {
        input_args[[i]][[j]] <- list(input_args[[i]][[j]])
      }
    }
  }
  if (any(sapply(input_args, \(x) !is.list(x))))
    stop("All inputs need to be lists of lists or vectors.")
  if (any(length(input_args[[1]]) != sapply(input_args, length)))
    stop("All input lists need to be of the same length.")

  scenarios <- list()
  scenario_idx <- list()
  for (i in seq_along(score_lists)){
    sublist <- lapply(input_args, \(x) x[[i]])
    names(sublist) <- argnames_wo_list
    scenarios[[i]] <- do.call(expand.grid, sublist)
    scenario_idx[i] <- nrow(scenarios[[i]])
  }
  scenarios <- do.call("rbind", scenarios)
  adestrOpts <- options()[startsWith(names(options()), "adestr_")]
  prog <- progressor(steps = nrow(scenarios))
  reslist <- future_apply(
    scenarios,
    MARGIN = 1L,
    FUN = \(x) {
      options(adestrOpts)
      res <- do.call(evaluate_estimator, unlist(x))
      prog(sprintf("Evaluating %s on %s for mu= %s.", toString(x$score), toString(x$estimator), format(x$mu)))
      res
    },
    future.seed = TRUE)
  res_split <- split(reslist, factor(rep(seq_along(scenario_idx), scenario_idx)))
  mergedlist <- list()
  for (reslist in res_split) {
    resnames <- unique(do.call("c", lapply(reslist, \(x)names(x@results))))
    splitlist <- lapply(resnames, \(x)list())
    names(splitlist) <- resnames
    for (i in seq_along(reslist)) {
      for (nm in resnames){
        tmpres <- reslist[[i]]
        if (nm %in% names(tmpres@results)) {
          row <- list(
            estimator = toString(tmpres@estimator),
            data_distribution = toString(tmpres@data_distribution),
            design = toString(tmpres@design),
            mu = tmpres@mu,
            sigma = tmpres@sigma,
            unnamed = tmpres@results[[nm]],
            error = tmpres@integrals[[nm]]$overall_integral$error,
            functionEvaluations = tmpres@integrals[[nm]]$overall_integral$functionEvaluations,
            idx = i
          )
          class(row) <- "data.frame"
          attr(row, "row.names") <- .set_row_names(length(row[[1]]))
          names(row)[[length(names(row))-3L]] <- nm
          names(row)[[length(names(row))-2L]] <- paste0("error_", nm)
          names(row)[[length(names(row))-1L]] <- paste0("functionEvaluations_", nm)
          splitlist[[nm]][[length(splitlist[[nm]])+1L]] <- row
        }
      }
    }
    for (nm in resnames){
      splitlist[[nm]] <- as.data.frame(do.call(rbind, splitlist[[nm]]))
    }
    merged <- splitlist[[1L]]
    for (i in seq_len(length(splitlist)-1L)){
      merged <- merge(merged, splitlist[[i+1L]], all=TRUE, by = c("estimator", "data_distribution", "design", "mu", "sigma"))
      merged$idx <- lapply(seq_along(merged$idx.x), \(i)unique(c(unlist(merged$idx.x[i]), unlist(merged$idx.y[i]))))
      merged$idx.x <- NULL
      merged$idx.y <- NULL
    }
    for (i in seq_len(length(merged) -1L)){
      merged[,i] <- unlist(merged[,i])
    }
    tmpidx <- merged$idx
    merged$idx <- NULL
    merged$EstimatorScoreResult <- lapply(tmpidx, \(idx) reslist[idx])
    class(merged$EstimatorScoreResult) <- c("EstimatorScoreResultList", class(merged$EstimatorScoreResult))
    mergedlist[[length(mergedlist)+1L]] <- merged
  }
  return(mergedlist)
}

