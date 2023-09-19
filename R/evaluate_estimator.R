#' Performance scores for point and interval estimators
#'
#' These classes encode various metrics which can be used to evaluate
#' the performance characteristics of point and interval estimators.
#'
#'
#' @slot label name of the performance score. Used in printing methods.
#'
#' @return an \code{EstimatorScore} object.
#' @export
#' @aliases EstimatorScore
#' @seealso \code{\link{evaluate_estimator}}
#'
#' @inherit evaluate_estimator examples
setClass("EstimatorScore", slots = c(label = "character"))
setClass("PointEstimatorScore", contains = "EstimatorScore")
setClass("IntervalEstimatorScore", contains = "EstimatorScore")
EstimatorScoreResult <- setClass("EstimatorScoreResult", slots = c(score = "list", estimator = "Estimator", data_distribution = "DataDistribution",
                                                                   design = "TwoStageDesign", mu = "ANY", sigma = "numeric",
                                                                   results = "list", integrals = "list"))
setClass("EstimatorScoreResultList", contains = "list")
EstimatorScoreResultList <- function(...) {
  r <- list(...)
  class(r) <- c("EstimatorScoreResultList", "list")
  r
}
setMethod("c", signature("EstimatorScoreResult"), definition =
            function(x, ...) {
              EstimatorScoreResultList(x, ...)
            })
setMethod("c", signature("EstimatorScoreResultList"), definition =
            function(x, ...) {
              EstimatorScoreResultList(x, ...)
            })

#' Evaluate performance characteristics of an estimator
#'
#' @param score performance measure to evaluate.
#' @param estimator object of class \code{PointEstimator}, \code{IntervalEstimator} or \code{PValue}.
#' @param data_distribution object of class \code{Normal} or \code{Student}.
#' @param use_full_twoarm_sampling_distribution logical indicating whether this estimator is intended to be used
#' with the full sampling distribution in a two-armed trial
#' @param design object of class \code{TwoStageDesign}.
#' @param true_parameter true value of the parameter (used e.g. when evaluating bias)
#' @param mu expected value of the underlying normal distribution
#' @param sigma assumed standard deviation.
#' @param tol relative tolerance
#' @param maxEval maximum number of iterations
#' @param absError absolute tolerance
#' @param exact logical indicating usage of exact n2 function.
#' @param early_futility_part include early futility part of integral
#' @param continuation_part include continuation part of integral
#' @param early_efficacy_part include early efficacy part of integral
#' @param conditional_integral treat integral as a conditional integral.
#'
#' @return \code{EstimatorScoreResult} object containing the performance characteristics of the estimator.
#' @export
#'
#' @examples
#' evaluate_estimator(
#'   score = MSE(),
#'   estimator = SampleMean(),
#'   data_distribution = Normal(FALSE),
#'   design = get_example_design(),
#'   mu = c(0.3),
#'   sigma = 1,
#'   exact = FALSE
#' )
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

#' @importFrom ggplot2 ggplot scale_x_continuous geom_line facet_wrap
#' @importFrom latex2exp TeX
setMethod("plot", signature = "EstimatorScoreResultList", definition =
            function(x, ...) {
              dat <- data.frame()
              for (estimatorscoreresult in x) {
                plot_list <- names(estimatorscoreresult@results)
                for (score in plot_list){
                  dat <- rbind(dat,
                               data.frame(mu = estimatorscoreresult@mu,
                                          Score = estimatorscoreresult@results[[score]],
                                          Estimator = toString(estimatorscoreresult@estimator),
                                          score_name = score
                                          )
                  )
                }
              }
              ggplot(data = dat, mapping = aes(x = mu, y = Score, col = Estimator)) +
                scale_x_continuous(name = TeX("$\\mu$")) +
                geom_line() +
                facet_wrap(vars(score_name), scales = "free_y")
            })
setMethod("plot", signature = "EstimatorScoreResult", definition =
            function(x, ...) {
              l <- EstimatorScoreResultList(x)
              plot(l, ...)
            })

#' @importFrom latex2exp TeX
plot_rejection_regions <- function(estimators, data_distribution, design, mu, sigma,  subdivisions = 100, ...){
  design <- TwoStageDesignWithCache(design)
  two_armed <- data_distribution@two_armed
  n1 <- n1(design, round = FALSE)
  se1 <- sigma_to_se(sigma, n1, two_armed)
  contl <- design@c1e - design@c1f
  minx <- design@c1f - 1.8*(1-2/(1+sqrt(5))) * contl
  maxx <- design@c1e + 2.2*(1-2/(1+sqrt(5))) * contl
  miny <- design@c1f - 4*(1-2/(1+sqrt(5))) * contl
  maxy <- design@c1e + (1-2/(1+sqrt(5))) * contl

  gridx <- seq(minx, maxx, length.out = subdivisions)
  gridy <- seq(miny, maxy, length.out = subdivisions)
  region <- expand.grid(T2 = gridy,
                        T1 = gridx)
  continuation_region <- region[design@c1f <= region$T1 & region$T1 <= design@c1e,]

  continuation_x_c2 <- seq(design@c1f, design@c1e, .01)
  alpha <- adoptr_alpha_shifted_design_kv(design, 0, 0, 0)

  c2_comb_list <- list()
  draw_line_3 <- FALSE
  for (i in seq_along(estimators)) {
    estimator <- estimators[[i]]
    stagewise_estimators <- get_stagewise_estimators(estimator = estimator,
                                                     data_distribution =  data_distribution,
                                                     use_full_twoarm_sampling_distribution = FALSE,
                                                     design = design, sigma = sigma, exact = FALSE)
    p1 <- stagewise_estimators[[1L]]
    p2 <- stagewise_estimators[[2L]]
    if (is(estimator, "StagewiseCombinationFunctionOrderingPValue")){
      continuation_y_c2 <- c2_extrapol(design, continuation_x_c2)
    } else{
      continuation_y_c2 <- implied_c2(design, continuation_x_c2, p2, sigma, two_armed, alpha)
    }
    x_line_3 <- NA
    if (p1(design = design, smean1 = z_to_smean(design@c1f, n1, sigma, two_armed), n1 = n1, sigma = sigma, two_armed = two_armed) > alpha) {
      yend1 <- maxy
    } else{
      yend1 <- miny
      draw_line_3 <- TRUE
      x_line_3 <- uniroot(\(x){
        p1(design = design, smean1 = z_to_smean(x, n1, sigma, two_armed), n1 = n1, sigma = sigma, two_armed = two_armed) - alpha
      },
      lower = -4, upper = 4, extendInt="yes")$root
    }
    if (p1(design = design, smean1 = z_to_smean(design@c1e, n1, sigma, two_armed), n1 = n1, sigma = sigma, two_armed = two_armed) < alpha) {
      yend2 <- miny
    } else{
      yend2 <- maxy
      draw_line_3 <- TRUE
      x_line_3<- uniroot(\(x){
        p1(design = design, smean1 = z_to_smean(x, n1, sigma, two_armed), n1 = n1, sigma = sigma, two_armed = two_armed) - alpha
      },
      lower = -4, upper = 4, extendInt="yes")$root
    }
    continuation_x_c2[c(1, length(continuation_x_c2))] <- continuation_x_c2[c(1, length(continuation_x_c2))]
    c2_comb <- data.frame(x = c(continuation_x_c2[1L],  continuation_x_c2, continuation_x_c2[length(continuation_x_c2)]),
                          y = c(yend1, continuation_y_c2, yend2),
                          x_line_3 = x_line_3,
                          Ordering = if (is(estimator, "NeymanPearsonOrderingPValue")) "Neyman-Pearson test ordering" else  substr(toTeX(estimator), 1, nchar(toTeX(estimator)) -8))
    c2_comb_list[[i]] <- c2_comb

  }
  c2_comb <- do.call("rbind", c2_comb_list)
  c2_comb$Ordering <- factor(c2_comb$Ordering, levels = unique(c2_comb$Ordering))

  limitsx <- c(min(gridx - .3), max(gridx + .3))
  limitsy <- c(min(gridy - .3), max(gridy + .3))
  plt <- ggplot() +
    geom_line(data = c2_comb, aes(x = .data$x, y = .data$y, color = .data$Ordering), linewidth=1) +
    scale_color_discrete(labels = sapply(as.character(unique(c2_comb$Ordering)), TeX)) +
    scale_x_continuous(name =TeX("$z_1$"),
                       limits = limitsx,
                       breaks = unique(round(seq(limitsx[1], limitsx[2], .5) )) ) +
    scale_y_continuous(name =TeX("$z_2$"),
                       limits = limitsy,
                       breaks = unique(round(seq(limitsy[1], limitsy[2], .5) ))) +
    theme_pubclean() +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("All rejection boundaries combined")
  if (draw_line_3) {
    c2_comb2 <- c2_comb[!is.na(c2_comb$x_line_3),]
    plt <- plt + geom_segment(data = c2_comb2, mapping = aes(x=x_line_3, xend=x_line_3, y=miny, yend = maxy, color=.data$Ordering), linewidth=1)
  }
  plt
}

setGeneric("plot_p", \(estimator, data_distribution, design, mu = 0, sigma, boundary_color="lightgreen", subdivisions = 100, ...) standardGeneric("plot_p"))
#' @import ggplot2 ggpubr latex2exp
setMethod("plot_p", signature("PValue"),
          function(estimator, data_distribution, design, mu, sigma, boundary_color, subdivisions, ...){
            design <- TwoStageDesignWithCache(design)
            two_armed <- data_distribution@two_armed
            n1 <- n1(design, round = FALSE)
            se1 <- sigma_to_se(sigma, n1, two_armed)
            contl <- design@c1e - design@c1f
            minx <- design@c1f - 1.8*(1-2/(1+sqrt(5))) * contl
            maxx <- design@c1e + 2.2*(1-2/(1+sqrt(5))) * contl
            miny <- design@c1f - 4*(1-2/(1+sqrt(5))) * contl
            maxy <- design@c1e + (1-2/(1+sqrt(5))) * contl

            gridx <- seq(minx, maxx, length.out = subdivisions)
            gridy <- seq(miny, maxy, length.out = subdivisions)
            region <- expand.grid(T2 = gridy,
                                  T1 = gridx)

            stagewise_estimators <- get_stagewise_estimators(estimator = estimator,
                                                             data_distribution =  data_distribution,
                                                             use_full_twoarm_sampling_distribution = FALSE,
                                                             design = design, sigma = sigma, exact = FALSE)
            p0 <- stagewise_estimators[[1L]]
            p1 <- stagewise_estimators[[2L]]

            p_futility <- sapply(gridx[gridx < design@c1f], \(T1) p0(design = design, smean1 = T1 * se1, n1 = n1, sigma = sigma, two_armed = two_armed))
            p_efficacy<- sapply(gridx[gridx > design@c1e], \(T1) p0(design = design, smean1 = T1 * se1, n1 = n1, sigma = sigma, two_armed = two_armed))

            continuation_region <- region[design@c1f <= region$T1 & region$T1 <= design@c1e,]
            p_continuation <- mapply(\(T1, T2) {
              n2 <- n2_extrapol(design, T1)
              se2 <- sigma / sqrt(n2)
              ret <- p1(
                design = design,
                smean1 = T1 * se1,
                smean2 = T2 * se2,
                n1 = n1,
                n2 = n2,
                sigma = sigma,
                two_armed = two_armed
              )
            },
            T1 = continuation_region$T1,
            T2 = continuation_region$T2)

            continuation_x_c2 <- seq(design@c1f, design@c1e, .01)
            alpha <- adoptr_alpha_shifted_design_kv(design, 0, 0, 0)
            if (is(estimator, "StagewiseOrderingPValue")){
              continuation_y_c2 <- c2_extrapol(design, continuation_x_c2)
            } else{
              continuation_y_c2 <- sapply(continuation_x_c2,
                                          \(x) uniroot(\(y) {
                                            n2 <- n2_extrapol(design, x)
                                            p1(design = design,
                                               smean1 = z_to_smean(x, n1, sigma, two_armed),
                                               smean2 = z_to_smean(y, n2, sigma, two_armed),
                                               n1 = n1, n2 = n2, sigma = sigma, two_armed = two_armed) - alpha
                                          } , lower = -4, upper = 4, extendInt = "yes")$root)
            }
            draw_line_3 <- FALSE
            if (p0(design = design, smean1 = z_to_smean(design@c1f, n1, sigma, two_armed), n1 = n1, sigma = sigma, two_armed = two_armed) > alpha) {
              yend1 <- maxy
            } else{
              yend1 <- miny
              draw_line_3 <- TRUE
              x_line_3 <- uniroot(\(x){
                p0(design = design, smean1 = z_to_smean(x, n1, sigma, two_armed), n1 = n1, sigma = sigma, two_armed = two_armed) - alpha
              },
              lower = -4, upper = 4, extendInt="yes")$root
            }
            if (p0(design = design, smean1 = z_to_smean(design@c1e, n1, sigma, two_armed), n1 = n1, sigma = sigma, two_armed = two_armed) < alpha) {
              yend2 <- miny
            } else{
              yend2 <- maxy
              draw_line_3 <- TRUE
              x_line_3<- uniroot(\(x){
                p0(design = design, smean1 = z_to_smean(x, n1, sigma, two_armed), n1 = n1, sigma = sigma, two_armed = two_armed) - alpha
              },
              lower = -4, upper = 4, extendInt="yes")$root
            }
            continuation_x_c2[c(1, length(continuation_x_c2))] <- continuation_x_c2[c(1, length(continuation_x_c2))]
            c2_comb <- data.frame(x = c(continuation_x_c2[1L],  continuation_x_c2, continuation_x_c2[length(continuation_x_c2)]),
                                  y = c(yend1, continuation_y_c2, yend2),
                                  PValue = estimator@label)
            p_comb <- rbind(
              cbind(
                T1 = gridx[gridx < design@c1f],
                T2 = rep(0, length(gridx[gridx < design@c1f])),
                p = p_futility
              ),
              cbind(continuation_region, p = p_continuation),
              cbind(
                T1 = gridx[gridx > design@c1e],
                T2 = rep(0, length(gridx[gridx > design@c1e])),
                p = p_efficacy
              )
            )

            p_comb <- cbind(region,
                            p = c(rep(p_futility, each = length(gridx)),
                                  p_continuation,
                                  rep(p_efficacy, each = length(gridx))
                            ))
            limitsx <- c(min(gridx - .3), max(gridx + .3))
            limitsy <- c(min(gridy - .3), max(gridy + .3))
            if (is(data_distribution, "Normal")) {
              xtxt <- TeX("$z_1$")
              ytxt <- TeX("$z_2$")
            } else if (is(data_distribution, "Student")) {
              xtxt <- TeX("$t_1$")
              ytxt <- TeX("$t_2$")
            } else {
             stop("unsupported data distribution")
            }
            plt <- ggplot() +
              # geom_tile(data = p_comb[p_comb$T1 < design@c1f & p_comb$T2==min(p_comb$T2),],
              #             aes(x = .data$`T1`, y = .data$`T2`, fill = .data$`p`, height=.25)) +
              # geom_tile(data = p_comb[p_comb$T1 >= design@c1f & p_comb$T1 <= design@c1e,],
              #           aes(x = .data$`T1`, y = .data$`T2`, fill = .data$`p`, height=.25)) +
              # geom_tile(data = p_comb[p_comb$T1 > design@c1e & p_comb$T2==max(p_comb$T2),],
              #             aes(x = .data$`T1`, y = .data$`T2`, fill = .data$`p`, height=.25)) +
              geom_tile(data = p_comb, aes(x = .data$`T1`, y = .data$`T2`, fill = .data$`p`)) +
              # geom_segment(data = c2_comb, mapping = aes(x=design@c1f, xend=design@c1f, y=.data[["y"]][[1]], yend = yend1), color="lightgreen", linewidth=1) +
              # geom_segment(data = c2_comb, mapping = aes(x=design@c1e, xend=design@c1e, y=.data[["y"]][[length(c2_comb[,1])]], yend = yend2), color="lightgreen", linewidth=1) +
              geom_line(data = c2_comb, aes(x = .data$x, y = .data$y), color = boundary_color, linewidth=1) +
              scale_color_manual() +
              scale_fill_gradient(low="blue", high="red") +
              scale_x_continuous(name =xtxt,
                                 limits = limitsx,
                                 breaks = unique(round(seq(limitsx[1], limitsx[2], .5) )) ) +
              scale_y_continuous(name =ytxt,
                                 limits = limitsy,
                                 breaks = unique(round(seq(limitsy[1], limitsy[2], .5) ))) +
              theme_pubclean() +
              theme(plot.title = element_text(hjust = 0.5), legend.key.width = unit(1, 'cm')) +
              ggtitle(TeX(toTeX(estimator)))
            if (draw_line_3) {
              plt <- plt + geom_segment(data = c2_comb, mapping = aes(x=x_line_3, xend=x_line_3, y=miny, yend = maxy), color = boundary_color, linewidth=1)
            }
            plt
          })


setGeneric("plot_sample_mean", \(data_distribution, design, mu, sigma, combine_components = TRUE, plot_treatment_group_if_twoarm = FALSE,
                                 p_limit = .0001, subdivisions = 100L,
                                 exact = FALSE, facets = 3L, ...) standardGeneric("plot_sample_mean"))
#' @import ggplot2
#' @importFrom forcats as_factor
setMethod("plot_sample_mean", signature("DataDistribution", "TwoStageDesign"),
          \(data_distribution, design, mu, sigma, combine_components, p_limit, subdivisions, exact, facets, ...) {
            two_armed <- data_distribution@two_armed
            n1 <- n1(design, round = FALSE)
            se1 <- sigma_to_se(sigma, n1, two_armed)
            smean_list <- list()
            for (m in mu){
              if (is(data_distribution, "Student")) {
                gridx <- suppressWarnings(seq(qt(1-p_limit, df = n1, ncp = m/se1, lower.tail = FALSE),
                                              qt(p_limit, df = n1, ncp = m/se1, lower.tail = FALSE),
                                              length.out = subdivisions)) * se1
                gridx <- gridx[!(abs(gridx)<.02)]
              } else {
                gridx <- seq(qnorm(1-p_limit, mean = m, sd = se1, lower.tail = FALSE),
                             qnorm(p_limit, mean = m, sd = se1, lower.tail = FALSE),
                             length.out = subdivisions)
              }
              if (plot_treatment_group_if_twoarm){
                y <- dsmeanT(data_distribution, design, smeanT = gridx, mu = m, sigma = sigma,
                             combine_components = combine_components, exact = exact)
              } else{
                y <- dsmean(data_distribution, design, smean = gridx, mu = m, sigma = sigma, two_armed = two_armed,
                            combine_components = combine_components, exact = exact)
              }

              if (combine_components)
                smean_list[[length(smean_list)+1L]] <- data.frame(smean = gridx, Density = y, mu = (paste0("mu == ",format(round(m, digits = 2)))))
            }
            if (combine_components){
              smean_dat <- do.call("rbind", smean_list)
              smean_dat$mu <- as_factor(smean_dat$mu)
              plt <- ggplot(data = smean_dat, aes(x = .data$`smean`, y = .data$`Density`)) +
                geom_line(size = 1) +
                scale_x_continuous(name = "Overall sample mean") +
                theme_pubclean()
              if (length(mu) > 1L) {
                plt <- plt + facet_wrap(vars(mu), labeller = label_parsed)
              }
              return(plt)
            } else {
              if (exact) {
                y[["0"]] <- y[["futility"]] + y[["efficacy"]]
                y$efficacy <- NULL
                y$futility <- NULL
                ys <- c(y["0"], y[seq(1L, length(y)-1L, length.out = facets)])
              }
              else
                ys <- c(y["futility"], y["continuation"], y["efficacy"])
              smean_dat <- data.frame(smean = rep(gridx, times = length(ys)),
                                      Density = do.call(c, ys),
                                      n = if (exact) paste0("n2 = ",rep(names(ys), each = length(gridx))) else rep(c("futility", "continuation", "efficacy"), each = length(gridx)) )
              smean_dat$n <- as_factor(smean_dat$n)
              plt <- ggplot(data = smean_dat, aes(x = .data$`smean`, y = .data$`Density`)) +
                geom_line(size = 1) +
                scale_x_continuous(name = "Sample mean") +
                facet_wrap(vars(n))
              return(plt)
            }
          })


#' Evaluate different scenarios in parallel
#'
#' This function takes a list of lists of scores, a list of lists of estimators,
#' and lists lists of various other design parameters. Each possible combination
#' of the elements of the respective sublists is then used to create a separate
#' scenario. These scenarios are than evaluated independelty of each other,
#' allowing for parallelization via the \code{future} framework.
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
#' @return list of data.frames containing the results for the respective scenarios.
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

