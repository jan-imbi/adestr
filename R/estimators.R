setClass("Statistic", slots = c(label = "character"))
setClass("Estimator", contains = "Statistic")

#' Point estimators
#'
#' @param g1 functional representation of the estimator in the early futility and efficacy regions.
#' @param g2 functional representation of the estimator in the continuation region.
#' @param label name of the estimator. Used in printing methods.
#'
#' @return An object of class \code{PointEstimator}.
#'
#' @export
#'
#' @examples
#' PointEstimator(g1 = \(smean1, ...) smean1,g2 = \(smean2, ...) smean2, label="My custom estimator")
setClass("PointEstimator", slots = c(g1 = "function", g2 = "function"), contains = "Estimator")
#' @rdname PointEstimator-class
#' @export
PointEstimator <- function(g1, g2, label) new("PointEstimator", g1 = g1, g2 = g2, label = label)
setClass("VirtualPointEstimator", contains = "PointEstimator")
#' @rdname PointEstimator-class
VirtualPointEstimator <- function() stop("Cannot create instance of class VirtualPointEstimator.")


#' P-values
#'
#' @param g1 functional representation of the p-value in the early futility and efficacy regions.
#' @param g2 functional representation of the p-value in the continuation region.
#' @param label name of the p-value. Used in printing methods.
#'
#' @return An object of class \code{PValue}.
#'
#' @export
#'
#' @examples
#' PValue(
#'   g1 = \(smean1, ...) runif(length(smean1)),
#'   g2 = \(smean2, ...) runif(length(smean2)),
#'   label="My custom p-value")
setClass("PValue", slots = c(g1 = "function", g2 = "function"),  contains = "Statistic")
#' @rdname PValue-class
#' @export
PValue <- function(g1, g2, label) new("PValue", g1 = g1, g2 = g2, label = label)
setClass("VirtualPValue", contains = "PValue")
#' @rdname PValue-class
VirtualPValue <- function() stop("Cannot create instance of class VirtualPValue.")


#' Interval estimators
#'
#' @param l1 functional representation of the lower boundary of the interval in the early futility and efficacy regions.
#' @param u1 functional representation of the upper boundary of the interval in the early futility and efficacy regions.
#' @param l2 functional representation of the lower boundary of the interval in the continuation region.
#' @param u2 functional representation of the upper boundary of the interval in the continuation region.
#' @param two_sided logical indicating whether the confidence interval is two-sided.
#' @param label name of the estimator. Used in printing methods.
#'
#' @return An object of class \code{IntervalEstimator}.
#'
#' @export
#'
#' @examples
#' IntervalEstimator(
#'   two_sided = FALSE,
#'   l1 = \(smean1, ...) smean1 - 1,
#'   u1 = \(smean1, ...) smean1 + 1,
#'   l2 = \(smean2, ...) smean2 - 1,
#'   u2 = \(smean2, ...) smean2 + 1,
#'   label="My custom p-value")
setClass("IntervalEstimator", slots = c(two_sided = "logical", l1 = "function", u1 = "function", l2 = "function", u2 = "function"), contains = "Estimator")
#' @rdname IntervalEstimator-class
#' @export
IntervalEstimator <- function(two_sided, l1, u1, l2, u2, label) new("IntervalEstimator", two_sided = two_sided, l1 = l1, u1 = u1, l2 = l2, u2 = u2, label = label)
setClass("VirtualIntervalEstimator", contains = "IntervalEstimator")
#' @rdname IntervalEstimator-class
#' @export
VirtualIntervalEstimator <- function() stop("Cannot create instance of class VirtualIntervalEstimator")



#' Conditional representations of an estimator or p-value
#'
#' This generic determines the functional representations of
#' point and interval estimators and p-values. The functions are
#' returned in two parts, one part to calculate the values conditional
#' on early futility or efficacy stops
#' (i.e. where no second stage mean and sample size is available),
#' and one conditional on continuation to the second stage.
#'
#' @inheritParams evaluate_estimator
#'
#' @return a list with the conditional functional representations
#' (one for each stage where the trial might end) of the estimator or p-value.
#' @export
#'
#' @examples
#' get_stagewise_estimators(
#'   estimator = SampleMean(),
#'   data_distribution = Normal(FALSE),
#'   use_full_twoarm_sampling_distribution = FALSE,
#'   design = get_example_design(),
#'   sigma = 1,
#'   exact = FALSE
#' )
setGeneric("get_stagewise_estimators", function(estimator,
                                                data_distribution,
                                                use_full_twoarm_sampling_distribution = FALSE,
                                                design,
                                                sigma,
                                                exact = FALSE) standardGeneric("get_stagewise_estimators"))
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("VirtualPointEstimator", "ANY"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            stop(sprintf("get_stagewise_estimators not implemented for estimator of class %s and data_distribution of class %s.", is(estimator)[[1L]],
                         is(data_distribution)[[1L]]))
          })
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("VirtualPValue", "ANY"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            stop(sprintf("get_stagewise_estimators not implemented for estimator of class %s and data_distribution of class %s.", is(estimator)[[1L]],
                         is(data_distribution)[[1L]]))
          })
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("VirtualIntervalEstimator", "ANY"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            stop(sprintf("get_stagewise_estimators not implemented for estimator of class %s and data_distribution of class %s.", is(estimator)[[1L]],
                         is(data_distribution)[[1L]]))
          })
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("PointEstimator", "Student"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            if (!existsMethod("get_stagewise_estimators", c(class(estimator), "Normal"))) {
              return(selectMethod("get_stagewise_estimators", signature=signature("PointEstimator", "DataDistribution"))(estimator,
                                                                                                                         data_distribution,
                                                                                                                         use_full_twoarm_sampling_distribution,
                                                                                                                         design,
                                                                                                                         sigma,
                                                                                                                         exact))
            }
            two_armed <- data_distribution@two_armed
            ests <- get_stagewise_estimators(estimator = estimator,
                                             data_distribution =  Normal(data_distribution@two_armed),
                                             use_full_twoarm_sampling_distribution = use_full_twoarm_sampling_distribution,
                                             design = design, sigma = sigma, exact = exact)
            g1_norm <- ests$g1
            g2_norm <- ests$g2
            g1 <- \(design, svar1, sigma, ...) mapply(g1_norm, sigma = sqrt(svar1), ..., MoreArgs = list(design = design))
            if (!two_armed) {
              g2_onearm <- \(design, smean1, svar1, smean2, svar2, n1, n2, sigma, two_armed, ...) mapply(g2_norm, smean1 = smean1, smean2 = smean2, n1 = n1, n2 = n2, sigma = sqrt(get_overall_svar_onearm(smean1, svar1, smean2, svar2, n1, n2)), two_armed = two_armed, ..., MoreArgs = list(design = design))
              return(list(g1 = g1,
                          g2 = g2_onearm))
            } else if (two_armed & use_full_twoarm_sampling_distribution) {
              g2_twoarm_full <- \(design, smean1, smean1T, svar1, smean2, smean2T, svar2, n1, n2, sigma,  ...) mapply(g2_norm, smean1 = smean1, smean1T = smean1T, smean2 = smean2, smean2T = smean2T, n1 = n1, n2 = n2, sigma = sqrt(get_overall_svar_twoarm(smean1, smean1T, svar1, smean2, smean2T, svar2, n1, n2)), ..., MoreArgs = list(design = design))
              return(list(g1 = g1,
                          g2 = g2_twoarm_full))
            } else {
              g2_twoarm <- \(design, svar1, svar2, n1, n2, sigma, two_armed, ...) mapply(g2_norm, n1 = n1, n2 = n2, sigma = sqrt(pool_svar1_svar2(svar1, svar2, n1, n2, two_armed)), two_armed = two_armed, ..., MoreArgs = list(design = design))
              return(list(g1 = g1,
                          g2 = g2_twoarm))
            }
            stop("Could not determine g2.")
          })
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("PValue", "Student"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            if (!existsMethod("get_stagewise_estimators", c(class(estimator), "Normal"))) {
              return(selectMethod("get_stagewise_estimators", signature=signature("PValue", "DataDistribution"))(estimator,
                                                                                                                 data_distribution,
                                                                                                                 use_full_twoarm_sampling_distribution,
                                                                                                                 design,
                                                                                                                 sigma,
                                                                                                                 exact))
            }
            two_armed <- data_distribution@two_armed
            ests <- get_stagewise_estimators(estimator = estimator,
                                             data_distribution =  Normal(data_distribution@two_armed),
                                             use_full_twoarm_sampling_distribution = use_full_twoarm_sampling_distribution,
                                             design = design, sigma = sigma, exact = exact)
            g1_norm <- ests$g1
            g2_norm <- ests$g2
            g1 <- \(design, svar1, sigma, ...) mapply(g1_norm, sigma = sqrt(svar1), ..., MoreArgs = list(design = design))
            if (!two_armed) {
              g2_onearm <- \(design, smean1, svar1, smean2, svar2, n1, n2, sigma, two_armed, ...) mapply(g2_norm, smean1 = smean1, smean2 = smean2, n1 = n1, n2 = n2, sigma = sqrt(get_overall_svar_onearm(smean1, svar1, smean2, svar2, n1, n2)), two_armed = two_armed, ..., MoreArgs = list(design = design))
              return(list(g1 = g1,
                          g2 = g2_onearm))
            } else if (two_armed & use_full_twoarm_sampling_distribution) {
              g2_twoarm_full <- \(design, smean1, smean1T, svar1, smean2, smean2T, svar2, n1, n2, sigma,  ...) mapply(g2_norm, smean1 = smean1, smean1T = smean1T, smean2 = smean2, smean2T = smean2T, n1 = n1, n2 = n2, sigma = sqrt(get_overall_svar_twoarm(smean1, smean1T, svar1, smean2, smean2T, svar2, n1, n2)), ..., MoreArgs = list(design = design))
              return(list(g1 = g1,
                          g2 = g2_twoarm_full))
            } else {
              g2_twoarm <- \(design, svar1, svar2, n1, n2, sigma, two_armed, ...) mapply(g2_norm, n1 = n1, n2 = n2, sigma = sqrt(pool_svar1_svar2(svar1, svar2, n1, n2, two_armed)), two_armed = two_armed, ..., MoreArgs = list(design = design))
              return(list(g1 = g1,
                          g2 = g2_twoarm))
            }
            stop("Could not determine g2.")
          })
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("IntervalEstimator", "Student"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            if (!existsMethod("get_stagewise_estimators", c(class(estimator), "Normal"))) {
              return(selectMethod("get_stagewise_estimators", signature=signature("IntervalEstimator", "DataDistribution"))(estimator,
                                                                                                                            data_distribution,
                                                                                                                            use_full_twoarm_sampling_distribution,
                                                                                                                            design,
                                                                                                                            sigma,
                                                                                                                            exact))
            }
            two_armed <- data_distribution@two_armed
            ests <- get_stagewise_estimators(estimator = estimator,
                                             data_distribution =  Normal(data_distribution@two_armed),
                                             use_full_twoarm_sampling_distribution = use_full_twoarm_sampling_distribution,
                                             design = design, sigma = sigma, exact = exact)
            l1_norm <- ests$l1
            l2_norm <- ests$l2
            u1_norm <- ests$u1
            u2_norm <- ests$u2
            l1 <- \(design, svar1, sigma, ...) mapply(l1_norm, sigma = sqrt(svar1), ..., MoreArgs = list(design = design))
            u1 <- \(design, svar1, sigma, ...) mapply(u1_norm, sigma = sqrt(svar1), ..., MoreArgs = list(design = design))
            if (!two_armed) {
              l2_onearm <- \(design, smean1, svar1, smean2, svar2, n1, n2, sigma, two_armed, ...) mapply(l2_norm, smean1 = smean1, smean2 = smean2, n1 = n1, n2 = n2, sigma = sqrt(get_overall_svar_onearm(smean1, svar1, smean2, svar2, n1, n2)), two_armed = two_armed, ..., MoreArgs = list(design = design))
              u2_onearm <- \(design, smean1, svar1, smean2, svar2, n1, n2, sigma, two_armed, ...) mapply(u2_norm, smean1 = smean1, smean2 = smean2, n1 = n1, n2 = n2, sigma = sqrt(get_overall_svar_onearm(smean1, svar1, smean2, svar2, n1, n2)), two_armed = two_armed, ..., MoreArgs = list(design = design))
              return(list(
                l1 = l1,
                u1 = u1,
                l2 = l2_onearm,
                u2 = u2_onearm
              ))
            } else if (two_armed & use_full_twoarm_sampling_distribution) {
              l2_twoarm_full <- \(design, smean1, smean1T, svar1, smean2, smean2T, svar2, n1, n2, sigma,  ...) mapply(l2_norm, smean1 = smean1, smean1T = smean1T, smean2 = smean2, smean2T = smean2T, n1 = n1, n2 = n2, sigma = sqrt(get_overall_svar_twoarm(smean1, smean1T, svar1, smean2, smean2T, svar2, n1, n2)), two_armed = two_armed, ..., MoreArgs = list(design=design))
              u2_twoarm_full <- \(design, smean1, smean1T, svar1, smean2, smean2T, svar2, n1, n2, sigma,  ...) mapply(u2_norm, smean1 = smean1, smean1T = smean1T, smean2 = smean2, smean2T = smean2T, n1 = n1, n2 = n2, sigma = sqrt(get_overall_svar_twoarm(smean1, smean1T, svar1, smean2, smean2T, svar2, n1, n2)), two_armed = two_armed, ..., MoreArgs = list(design=design))
              return(list(
                l1 = l1,
                u1 = u1,
                l2 = l2_twoarm_full,
                u2 = u2_twoarm_full
              ))
            } else {
              l2_twoarm <- \(design, svar1, svar2, n1, n2, sigma, two_armed, ...) mapply(l2_norm, n1 = n1, n2 = n2, sigma = sqrt(pool_svar1_svar2(svar1, svar2, n1, n2, two_armed)), two_armed = two_armed, ..., MoreArgs = list(design=design))
              u2_twoarm <- \(design, svar1, svar2, n1, n2, sigma, two_armed, ...) mapply(u2_norm, n1 = n1, n2 = n2, sigma = sqrt(pool_svar1_svar2(svar1, svar2, n1, n2, two_armed)), two_armed = two_armed, ..., MoreArgs = list(design=design))
              return(list(
                l1 = l1,
                u1 = u1,
                l2 = l2_twoarm,
                u2 = u2_twoarm
              ))
            }
            stop("Could not determine g2.")
          })
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("VirtualPointEstimator", "Student"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {            if (!existsMethod("get_stagewise_estimators", c(class(estimator), "Normal"))) {
              return(selectMethod("get_stagewise_estimators", signature=signature("VirtualPointEstimator", "ANY"))(estimator,
                                                                                                                   data_distribution,
                                                                                                                   use_full_twoarm_sampling_distribution,
                                                                                                                   design,
                                                                                                                   sigma,
                                                                                                                   exact))
            } else{
              return(selectMethod("get_stagewise_estimators", signature=signature("PointEstimator", "Student"))(estimator,
                                                                                                                data_distribution,
                                                                                                                use_full_twoarm_sampling_distribution,
                                                                                                                design,
                                                                                                                sigma,
                                                                                                                exact))
            }
          })
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("VirtualIntervalEstimator", "Student"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            if (!existsMethod("get_stagewise_estimators", c(class(estimator), "Normal"))) {
              return(selectMethod("get_stagewise_estimators", signature=signature("VirtualIntervalEstimator", "ANY"))(estimator,
                                                                                                                      data_distribution,
                                                                                                                      use_full_twoarm_sampling_distribution,
                                                                                                                      design,
                                                                                                                      sigma,
                                                                                                                      exact))
            } else{
              return(selectMethod("get_stagewise_estimators", signature=signature("IntervalEstimator", "Student"))(estimator,
                                                                                                                   data_distribution,
                                                                                                                   use_full_twoarm_sampling_distribution,
                                                                                                                   design,
                                                                                                                   sigma,
                                                                                                                   exact))
            }
          })
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("VirtualPValue", "Student"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            if (!existsMethod("get_stagewise_estimators", c(class(estimator), "Normal"))) {
              return(selectMethod("get_stagewise_estimators", signature=signature("VirtualPValue", "ANY"))(estimator,
                                                                                                           data_distribution,
                                                                                                           use_full_twoarm_sampling_distribution,
                                                                                                           design,
                                                                                                           sigma,
                                                                                                           exact))
            } else{
              return(selectMethod("get_stagewise_estimators", signature=signature("PValue", "Student"))(estimator,
                                                                                                        data_distribution,
                                                                                                        use_full_twoarm_sampling_distribution,
                                                                                                        design,
                                                                                                        sigma,
                                                                                                        exact))
            }
          })
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("PointEstimator", "DataDistribution"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            list(g1 = estimator@g1,
                 g2 = estimator@g2)
          })
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("PValue", "DataDistribution"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            list(g1 = estimator@g1,
                 g2 = estimator@g2)
          })
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("IntervalEstimator", "DataDistribution"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            list(l1 = estimator@l1,
                 u1 = estimator@u1,
                 l2 = estimator@l2,
                 u2 = estimator@u2)
          })


setClass("SampleMean", contains = "PointEstimator")
#' @rdname PointEstimator-class
#' @export
SampleMean <- function() new("PointEstimator", g1 = \(smean1, ...) smean1, g2 = \(smean1, smean2, n1, n2, ...) (n1 * smean1 + n2 * smean2) / (n1 + n2), label = "Sample mean")
setClass("FirstStageSampleMean", contains = "PointEstimator")
#' @rdname PointEstimator-class
#' @export
FirstStageSampleMean <- function() new("PointEstimator", g1 = \(smean1, ...) smean1, g2 = \(smean1, ...) smean1, label = "First-stage sample mean")
setClass("WeightedSampleMean", contains = "PointEstimator")

#' @rdname PointEstimator-class
#' @param w1 weight of the first-stage sample mean.
#' @importFrom scales percent
#' @export
WeightedSampleMean <- function(w1=.5) new("PointEstimator", g1 = \(smean1, ...) smean1, g2 = \(smean1, smean2, n1, n2, ...) w1 * smean1 + (1-w1) * smean2,
                                          label = paste0("WeightedSampleMean(w1=", percent(w1), ")", collapse=""))

setClass("AdaptivelyWeightedSampleMean", contains = "VirtualPointEstimator", slots = c(w1 = "numeric"))
#' @rdname PointEstimator-class
#' @importFrom scales percent
AdaptivelyWeightedSampleMean <- function(w1 = 1/sqrt(2)) {
  new(
    "AdaptivelyWeightedSampleMean",
    w1 = w1,
    label = paste0(
      "AdaptivelyWeightedSampleMean(w1^2=",
      scales::percent(w1^2),
      ")",
      collapse = ""
    )
  )
}
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("AdaptivelyWeightedSampleMean", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            w1 <- estimator@w1
            g1 = \(smean1, ...) smean1
            g2 = \(smean1, smean2, n1, n2, sigma, two_armed, ...) {
              se1 <- sigma_to_se(sigma, n1, two_armed)
              se2 <- sigma_to_se(sigma, n2, two_armed)
              tau <- (w1/se1) / (w1/se1 + sqrt(1-w1^2)/se2)
              tau * smean1 + (1 - tau) * smean2
            }
            list(g1 = g1, g2 = g2)
})

setClass("MinimizePeakVariance", contains = "VirtualPointEstimator")
#' @rdname PointEstimator-class
#' @export
MinimizePeakVariance <- function() new("MinimizePeakVariance", label = "AWSM minimizing peak variance")
#' @rdname get_stagewise_estimators
#' @importFrom stats optimize
setMethod("get_stagewise_estimators", signature("MinimizePeakVariance", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            get_ess <- function(mu) {
              H <- PointMassPrior(mu, 1)
              ess <- ExpectedSampleSize(data_distribution, H)
              evaluate(ess, design, optimization=TRUE)
            }
            max_ess_mu <- optimize(get_ess, z_to_smean(c(design@c1f, design@c1e), n1(design, round=FALSE), sigma, data_distribution@two_armed), maximum = TRUE)$maximum
            get_var <- function(w1) {
              est <- AdaptivelyWeightedSampleMean(w1)
              evaluate_estimator(score = Variance(),
                                 estimator = est,
                                 data_distribution = data_distribution,
                                 use_full_twoarm_sampling_distribution = use_full_twoarm_sampling_distribution,
                                 design = design,
                                 mu = max_ess_mu,
                                 sigma = sigma,
                                 tol = getOption("adestr_tol_inner", default = .adestr_options[["adestr_tol_inner"]]),
                                 maxEval = getOption("adestr_maxEval_inner", default = .adestr_options[["adestr_maxEval_inner"]]),
                                 absError = getOption("adestr_absError_inner", default = .adestr_options[["adestr_absError_inner"]]),
                                 exact = exact
                                 )@results$Variance
            }
            min_var_w1 <- optimize(get_var, c(1e-10, 1-1e-10))$minimum
            get_stagewise_estimators(estimator = AdaptivelyWeightedSampleMean(min_var_w1),
                                                             data_distribution =  data_distribution,
                                                             use_full_twoarm_sampling_distribution = use_full_twoarm_sampling_distribution,
                                                             design = design, sigma = sigma, exact = exact)
          })



setClass("BiasReduced", contains = "VirtualPointEstimator", slots = c(iterations = "numeric"))
#' @rdname PointEstimator-class
#' @param iterations number of bias reduction iterations. Defaults to 1.
#' @export
BiasReduced <- function(iterations = 1L) new("BiasReduced", iterations = iterations,
                                                         label = paste0("Bias reduced MLE (iterations=", iterations,")", collapse = ""))
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("BiasReduced", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            mle_expectation <- Vectorize(\(design, mu_plugin, sigma, two_armed, ...) {
              int_kv(design = design,
                     g1 = \(smean1, ...) smean1,
                     g2 = \(smean1, smean2, n1, n2, ...) (n1 * smean1 + n2 * smean2) / (n1 + n2),
                     mu = mu_plugin,
                     sigma = sigma,
                     two_armed = two_armed,
                     tol = getOption("adestr_tol_inner", default = .adestr_options[["adestr_tol_inner"]]),
                     maxEval = getOption("adestr_maxEval_inner", default = .adestr_options[["adestr_maxEval_inner"]]),
                     absError = getOption("adestr_absError_inner", default = .adestr_options[["adestr_absError_inner"]]),
                     exact = FALSE)$overall_integral$integral}, c("mu_plugin"))
            diff <- Vectorize(\(design, mu_plugin, sigma, two_armed, ...) {
              int_kv(design = design,
                     g1 = \(smean1, n1, sigma, two_armed, ...) {
                       se1sq <- sigma^2 * (1L + two_armed) / n1
                       (smean1) * ( (smean1 - mu_plugin) /se1sq)
                     },
                     g2 = \(smean1, smean2, n1, n2, sigma, two_armed, ...) {
                       se1sq <- sigma^2 * (1L + two_armed) / n1
                       se2sq <- sigma^2 * (1L + two_armed) / n2
                       ((n1 * smean1 + n2 * smean2) / (n1 + n2)) * ( (smean1 - mu_plugin) /se1sq + (smean2 - mu_plugin) /se2sq )
                     },
                     mu = mu_plugin,
                     sigma = sigma,
                     two_armed = two_armed,
                     tol = getOption("adestr_tol_inner", default = .adestr_options[["adestr_tol_inner"]]),
                     maxEval = getOption("adestr_maxEval_inner", default = .adestr_options[["adestr_maxEval_inner"]]),
                     absError = getOption("adestr_absError_inner", default = .adestr_options[["adestr_absError_inner"]]),
                     exact = FALSE)$overall_integral$integral}, c("mu_plugin"))
            g1 <- \(design, smean1, sigma, two_armed, ...) {
              estimate0 <- smean1
              estimate <- estimate0
              for (i in seq_len(estimator@iterations)){
                denom <- 1 + diff(design, mu_plugin = estimate, sigma, two_armed, ...)
                estimate <- estimate + ((estimate0 - estimate) -  (mle_expectation(design, mu_plugin = estimate, sigma, two_armed, ...) - estimate0) ) / denom
              }
              return(estimate)
            }
            g2 <- \(design, smean1, smean2, n1, n2, sigma, two_armed, ...) {
              estimate0 <- (smean1 * n1 + smean2 * n2) / (n1 + n2)
              estimate <- estimate0
              for (i in seq_len(estimator@iterations)){
                denom <- 1 + diff(design, mu_plugin = estimate, sigma, two_armed, ...)
                estimate <- estimate + ((estimate0 - estimate) -  (mle_expectation(design, mu_plugin = estimate, sigma, two_armed, ...) - estimate0)) / denom
              }
              return(estimate)
            }
            list(g1 = g1,
                 g2 = g2)
          })

rb1_kv <- function(smean1, n1, ...){
  smean1
}
rb2_kv <- function(smean1, smean2, n1, n2, mu, sigma, two_armed, preimage,
                  tol = getOption("adestr_tol_inner", default = .adestr_options[["adestr_tol_inner"]]),
                  maxEval = getOption("adestr_maxEval_inner", default = .adestr_options[["adestr_maxEval_inner"]]),
                  absError = getOption("adestr_absError_inner", default = .adestr_options[["adestr_absError_inner"]]), ...){
  if (missing(mu))
    mu <- 0
  smean <- (smean1*n1 + smean2 * n2) / (n1 + n2)
  denom <- .hcubature(
    \(x) {
      smean1_prime <- z_to_smean(x[1L,,drop=FALSE], n1, sigma, two_armed)
      z <- rbind(x[1L,,drop=FALSE], smean_to_z((smean * (n1 + n2) - smean1_prime * n1) / n2, n2, sigma, two_armed))
      mf2_kv(z, n1, n2, mu, sigma, two_armed)
    },
    lowerLimit = preimage[1L],
    upperLimit = preimage[2L],
    tol = tol,
    maxEval = maxEval,
    absError = absError,
    vectorInterface = TRUE)$integral
  .hcubature(
    \(x) {
      smean1_prime <- z_to_smean(x[1L,,drop=FALSE], n1, sigma, two_armed)
      z <- rbind(x[1L,,drop=FALSE], smean_to_z((smean * (n1 + n2) - smean1_prime * n1) / n2, n2, sigma, two_armed))
      smean1_prime*mf2_kv(z, n1, n2, mu, sigma, two_armed)
    },
    lowerLimit = preimage[1L],
    upperLimit = preimage[2L],
    tol = tol,
    maxEval = maxEval,
    absError = absError,
    vectorInterface = TRUE)$integral / denom
}

newrb2_kv <- function(smean, n1, n2, mu, sigma, two_armed, preimage,
                   tol = getOption("adestr_tol_inner", default = .adestr_options[["adestr_tol_inner"]]),
                   maxEval = getOption("adestr_maxEval_inner", default = .adestr_options[["adestr_maxEval_inner"]]),
                   absError = getOption("adestr_absError_inner", default = .adestr_options[["adestr_absError_inner"]])){
  denom <- .hcubature(
    \(x) {
      smean1_prime <- z_to_smean(x[1L,,drop=FALSE], n1, sigma, two_armed)
      z <- rbind(x[1L,,drop=FALSE], smean_to_z((smean * (n1 + n2) - smean1_prime * n1) / n2, n2, sigma, two_armed))
      mf2_kv(z, n1, n2, mu, sigma, two_armed)
    },
    lowerLimit = preimage[1L],
    upperLimit = preimage[2L],
    tol = tol,
    maxEval = maxEval,
    absError = absError,
    vectorInterface = TRUE)$integral
  .hcubature(
    \(x) {
      smean1_prime <- z_to_smean(x[1L,,drop=FALSE], n1, sigma, two_armed)
      z <- rbind(x[1L,,drop=FALSE], smean_to_z((smean * (n1 + n2) - smean1_prime * n1) / n2, n2, sigma, two_armed))
      smean1_prime*mf2_kv(z, n1, n2, mu, sigma, two_armed)
    },
    lowerLimit = preimage[1L],
    upperLimit = preimage[2L],
    tol = tol,
    maxEval = maxEval,
    absError = absError,
    vectorInterface = TRUE)$integral / denom
}
crb1_kv <- function(smean1, n1, ...){
  NA_real_
}
crb2_kv <- function(smean1, smean2, n1, n2, mu, sigma, two_armed, preimage,
                   tol = getOption("adestr_tol_inner", default = .adestr_options[["adestr_tol_inner"]]),
                   maxEval = getOption("adestr_maxEval_inner", default = .adestr_options[["adestr_maxEval_inner"]]),
                   absError = getOption("adestr_absError_inner", default = .adestr_options[["adestr_absError_inner"]]),
                   ...){
  if (missing(mu))
    mu <- 0
  smean <- (smean1*n1 + smean2 * n2) / (n1 + n2)
  denom <- .hcubature(
    \(x) {
      smean1_prime <- z_to_smean(x[1L,,drop=FALSE], n1, sigma, two_armed)
      smean2_prime <- (smean * (n1 + n2) - smean1_prime * n1) / n2
      z <- rbind(x[1L,,drop=FALSE], smean_to_z(smean2_prime, n2, sigma, two_armed))
      mf2_kv(z, n1, n2, mu, sigma, two_armed)
    },
    lowerLimit = preimage[1L],
    upperLimit = preimage[2L],
    tol = tol,
    maxEval = maxEval,
    absError = absError,
    vectorInterface = TRUE)$integral
  .hcubature(
    \(x) {
      smean1_prime <- z_to_smean(x[1L,,drop=FALSE], n1, sigma, two_armed)
      smean2_prime <- (smean * (n1 + n2) - smean1_prime * n1) / n2
      z <- rbind(x[1L,,drop=FALSE], smean_to_z(smean2_prime, n2, sigma, two_armed))
      smean2_prime*mf2_kv(z, n1, n2, mu, sigma, two_armed)
    },
    lowerLimit = preimage[1L],
    upperLimit = preimage[2L],
    tol = tol,
    maxEval = maxEval,
    absError = absError,
    vectorInterface = TRUE)$integral / denom
}

pseudorb2_kv <- function(design, smean1, smean2, n1, n2, mu, sigma, two_armed,
                   tol = getOption("adestr_tol_inner", default = .adestr_options[["adestr_tol_inner"]]),
                   maxEval = getOption("adestr_maxEval_inner", default = .adestr_options[["adestr_maxEval_inner"]]),
                   absError = getOption("adestr_absError_inner", default = .adestr_options[["adestr_absError_inner"]]),
                   ...){
  if (missing(mu))
    mu <- 0
  smean <- (smean1*n1 + smean2 * n2) / (n1 + n2)
  denom <- .hcubature(
    \(x) {
      smean1_prime <- z_to_smean(x[1L,,drop=FALSE], n1, sigma, two_armed)
      z <- rbind(x[1L,,drop=FALSE], smean_to_z((smean * (n1 + n2) - smean1_prime * n1) / n2, n2, sigma, two_armed))
      mf2_kv(z, n1, n2, mu, sigma, two_armed)
    },
    lowerLimit = design@c1f,
    upperLimit = design@c1e,
    tol = tol,
    maxEval = maxEval,
    absError = absError,
    vectorInterface = TRUE)$integral
  .hcubature(
    \(x) {
      smean1_prime <- z_to_smean(x[1L,,drop=FALSE], n1, sigma, two_armed)
      z <- rbind(x[1L,,drop=FALSE], smean_to_z((smean * (n1 + n2) - smean1_prime * n1) / n2, n2, sigma, two_armed))
      smean1_prime*mf2_kv(z, n1, n2, mu, sigma, two_armed)
    },
    lowerLimit = design@c1f,
    upperLimit = design@c1e,
    tol = tol,
    maxEval = maxEval,
    absError = absError,
    vectorInterface = TRUE)$integral / denom
}

setClass("RaoBlackwell", contains = "VirtualPointEstimator")
#' @rdname PointEstimator-class
#' @export
RaoBlackwell <- function() new("RaoBlackwell", label = "Rao-Blackwellized")
#' @rdname get_stagewise_estimators
#' @importFrom utils tail
setMethod("get_stagewise_estimators", signature("RaoBlackwell", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            if (exact) {
              return(list(g1 = rb1_kv,
                          g2 = Vectorize(rb2_kv, c("smean1", "smean2")) ))
            } else {
              sample_space <- n2_preimage(design = design, sigma = sigma, two_armed = data_distribution@two_armed)
              intervals <- c(sapply(sample_space, \(x)x$preimage[1]), tail(sample_space)$preimage[2])
              regions <- lapply(sample_space, \(x)x$preimage)
              ns <- sapply(sample_space, \(x)x$n)
              rb2_kv_vec <- Vectorize(rb2_kv, c("smean1", "smean2", "n2", "preimage"))
              g2 <- function(smean1, smean2, n1, n2, mu, sigma, two_armed, ...) {
                if (missing(mu))
                  mu <- 0
                indices <- findInterval(smean_to_z(smean1, n1, sigma, two_armed) , intervals)
                if (any(indices >length(sample_space)))
                  stop("Index for smean1 region out of bounds.")
                preimages <- regions[indices]
                #                           oder ns[indices] ?
                rb2_kv_vec(smean1, smean2, n1, n2, mu, sigma, two_armed, preimages)
              }
              return(list(g1 = rb1_kv,
                          g2 = g2 ))
            }
          })

setClass("PseudoRaoBlackwell", contains = "VirtualPointEstimator")
#' @rdname PointEstimator-class
#' @export
PseudoRaoBlackwell <- function() new("PseudoRaoBlackwell", label = "Pseudo Rao-Blackwellized")
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("PseudoRaoBlackwell", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            return(list(g1 = rb1_kv,
                        g2 = Vectorize(pseudorb2_kv, c("smean1", "smean2", "n2")) ))
          })

setClass("RepeatedCI", contains = "VirtualIntervalEstimator")
#' @rdname IntervalEstimator-class
#' @export
RepeatedCI <- function(two_sided = TRUE) new("RepeatedCI", two_sided = two_sided, label = "Repeated CI")
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("RepeatedCI", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            if (estimator@two_sided){
              list(l1 = Vectorize(rcil1_kv, c("smean1")),
                   u1 = Vectorize(rciu1_kv, c("smean1")),
                   l2 = Vectorize(rcil2_kv, c("smean1", "smean2", "n2")),
                   u2 = Vectorize(rciu2_kv, c("smean1", "smean2", "n2")))
            } else {
              list(l1 = Vectorize(rcil1_kv, c("smean1")),
                   u1 = \(...) Inf,
                   l2 = Vectorize(rcil2_kv, c("smean1", "smean2", "n2")),
                   u2 = \(...) Inf)
            }
          })

combf_kv <- function(design, smean1, smean2, n1, n2, mu, sigma, two_armed, ...) {
  se1 <- sigma_to_se(sigma, n1, two_armed)
  se2 <- sigma_to_se(sigma, n2, two_armed)
  y1 <- (smean1 - mu) / se1
  (smean2 - mu) / se2 - c2_extrapol(design, y1)
}
rcil1_kv <- function(design, smean1, n1, sigma, two_armed, ...){
  se1 <- sigma_to_se(sigma, n1, two_armed)
  smean1 - design@c1e * se1
}
rciu1_kv <- function(design, smean1, n1, sigma, two_armed, ...){
  se1 <- sigma_to_se(sigma, n1, two_armed)
  smean1 + design@c1e * se1
}
rcil2_kv <- function(design, smean1, smean2, n1, n2, sigma, two_armed, tol = getOption("adestr_tol_roots", default = .adestr_options[["adestr_tol_roots"]]), maxiter = getOption("adestr_maxiter_roots", default = .adestr_options[["adestr_maxiter_roots"]]), ...){
  se1 <- sigma_to_se(sigma, n1, two_armed)
  lb <- smean1 - design@c1e * se1
  ub <- smean1 - design@c1f * se1
  f.lb <- combf_kv(design, smean1, smean2, n1, n2, lb, sigma, two_armed)
  f.ub <- combf_kv(design, smean1, smean2, n1, n2, ub, sigma, two_armed)
  if (f.ub >= 0){
    root2 <- ub
  } else if (sign(f.ub)==sign(f.lb)) {
    root2 <- lb
  } else {
    root2 <- uniroot(
      \(x) combf_kv(design, smean1, smean2, n1, n2, x, sigma, two_armed),
      c(lb, ub),
      extendInt = "no",
      tol = tol,
      maxiter = maxiter
    )$root
  }
  root2
}
rciu2_kv <- function(design, smean1, smean2, n1, n2, sigma, two_armed, tol = getOption("adestr_tol_roots", default = .adestr_options[["adestr_tol_roots"]]), maxiter = getOption("adestr_maxiter_roots", default = .adestr_options[["adestr_maxiter_roots"]]), ...) -rcil2_kv(design = design, smean1 = -smean1, smean2 = -smean2, n1 = n1, n2 = n2, sigma = sigma, two_armed = two_armed, tol = tol, maxiter = maxiter, ...)


adoptr_alpha_shifted_design_kv <- function(design, shiftc1f, shiftc1e, shiftc2){
  design <- TwoStageDesignWithCache(design)
  pr_es1 <- pnorm(design@c1e + shiftc1e, mean = 0, sd = 1, lower.tail = FALSE)
  if (design@c1e + shiftc1e > design@c1f + shiftc1f + .Machine$double.eps){
    pr_es2 <-  .hcubature(f = \(x) matrix(pnorm(c2_extrapol(design, x[1L,]) + shiftc2, lower.tail = FALSE)*dnorm(x[1L,]), nrow=1L),
                         lowerLimit = design@c1f + shiftc1f,
                         upperLimit = design@c1e + shiftc1e,
                         tol = getOption("adestr_tol_inner", default = .adestr_options[["adestr_tol_inner"]]),
                         maxEval = getOption("adestr_maxEval_inner", default = .adestr_options[["adestr_maxEval_inner"]]),
                         absError = getOption("adestr_absError_inner", default = .adestr_options[["adestr_absError_inner"]]),
                         vectorInterface = TRUE)$integral
  } else{
    pr_es2 <- 0
  }
  pr_es1 + pr_es2
}

setClass("LinearShiftRepeatedPValue", slots =  c(wc1f="numeric", wc1e="numeric", wc2="numeric"), contains = "VirtualPValue")
#' @rdname PValue-class
#' @param wc1f slope of futility boundary change.
#' @param wc1e slope of efficacy boundary change.
#' @param wc2 slope of c2 boundary change.
#' @export
LinearShiftRepeatedPValue <- function(wc1f=0, wc1e=1/2, wc2=1/2) new("LinearShiftRepeatedPValue", wc1f=wc1f, wc1e=wc1e, wc2=wc2, label = "Repeated p value")
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("LinearShiftRepeatedPValue", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            list(g1 = Vectorize(\(design, smean1, n1, sigma, two_armed, ...) rp1_kv(design, smean1, sigma, wc1f = estimator@wc1f, wc1e= estimator@wc1e, wc2 = estimator@wc2, ...), c("smean1") ),
                 g2 = Vectorize(\(design, smean1, smean2, n1, n2, sigma, two_armed, ...) rp2_kv(design, smean1, sigma, wc1f=estimator@wc1f, wc1e=estimator@wc1e, wc2=estimator@wc2, ...), c("smean1", "smean2", "n2")))
          })
rp1_kv <- function(design, smean1, n1, sigma, two_armed, wc1f=0, wc1e=1/2, wc2=1/2, ...) {
  se1 <- sigma_to_se(sigma, n1, two_armed)
  diff <- design@c1e - smean1 / se1
  adoptr_alpha_shifted_design_kv(design = design,
                                shiftc1f = diff*wc1f,
                                shiftc1e = diff*wc1e,
                                shiftc2  = diff*wc2)
}
rp2_kv <- function(design, smean1, smean2, n1, n2, sigma, two_armed, wc1f=0, wc1e=1/2, wc2=1/2, ...) {
  diff <- -combf_kv(design, smean1, smean2, n1, n2, mu=0, sigma, two_armed)
  adoptr_alpha_shifted_design_kv(design = design,
                                shiftc1f = diff*wc1f,
                                shiftc1e = diff*wc1e,
                                shiftc2  = diff*wc2)
}

p_ml <- function(design, smean, n, mu, sigma, two_armed, tol = getOption("adestr_tol_inner", default = .adestr_options[["adestr_tol_inner"]]), maxEval = getOption("adestr_maxEval_inner", default = .adestr_options[["adestr_maxEval_inner"]]), absError = getOption("adestr_absError_inner", default = .adestr_options[["adestr_absError_inner"]]), ...) {
  design <- TwoStageDesignWithCache(design)
  n1 <- n1(design, round=FALSE)
  se1 <- sigma_to_se(sigma, n1, two_armed)
  muse1 <- mu/se1
  # smean1_sol <- smean
  z1 <- smean_to_z(smean, n1, sigma, two_armed)
  if (z1 < design@c1f) {
    p_part_first_stage <- pnorm(design@c1f, mean = muse1) - pnorm(z1, mean = muse1) + pnorm(design@c1e, mean = muse1, lower.tail = FALSE)
  } else if (z1 > design@c1e) {
    p_part_first_stage <- pnorm(z1, mean = muse1, lower.tail = FALSE)
  } else {
    p_part_first_stage <- pnorm(design@c1e, mean = muse1, lower.tail = FALSE)
  }
  p_part_second_stage <- .hcubature(f = \(x) {
    n2 <- n2_extrapol(design, x[1L,])
    se2 <- sigma * sqrt((1L + two_armed) /  n2)
    smean1_prime <- z_to_smean(x[1L,], n1, sigma, two_armed)
    f1_kv(x[1L,,drop=FALSE], n1, mu, sigma, two_armed) *
      pnorm( ((smean * (n1 + n2) - n1 * smean1_prime) / n2) /se2,
             mean = mu/se2,lower.tail = FALSE)
  },
  lowerLimit = c(design@c1f),
  upperLimit = c(design@c1e),
  tol = tol,
  maxEval = maxEval,
  absError = absError,
  vectorInterface = TRUE)$integral
  return(p_part_first_stage + p_part_second_stage)
}
p1_ml <- function(design, smean1, n1, mu, sigma, two_armed, ...) p_ml(design, smean1, n1, mu, sigma, two_armed, ...)
p2_ml <- function(design, smean1, smean2, n1, n2, mu, sigma, two_armed, ...) p_ml(design, smeans_to_smean(smean1, smean2, n1, n2), (n1 + n2), mu, sigma, two_armed, ...)

p_lr <- function(design, smean, n, mu, sigma, two_armed, tol = getOption("adestr_tol_inner", default = .adestr_options[["adestr_tol_inner"]]), maxEval = getOption("adestr_maxEval_inner", default = .adestr_options[["adestr_maxEval_inner"]]), absError = getOption("adestr_absError_inner", default = .adestr_options[["adestr_absError_inner"]]), ...) {
  design <- TwoStageDesignWithCache(design)
  n1 <- n1(design, round=FALSE)
  se1 <- sigma_to_se(sigma, n1, two_armed)
  muse1 <- mu/se1
  smean1_sol <- (smean * sqrt(n) + mu * sqrt(n1) - mu * sqrt(n)) / sqrt(n1)
  z1 <- smean_to_z(smean1_sol, n1, sigma, two_armed)
  if (z1 < design@c1f) {
    p_part_first_stage <- pnorm(design@c1f, mean = muse1) - pnorm(z1, mean = muse1) + pnorm(design@c1e, mean = muse1, lower.tail = FALSE)
  } else if (z1 > design@c1e) {
    p_part_first_stage <- pnorm(z1, mean = muse1, lower.tail = FALSE)
  } else {
    p_part_first_stage <- pnorm(design@c1e, mean = muse1, lower.tail = FALSE)
  }
  sqrtn <- sqrt(n)
  p_part_second_stage <- .hcubature(f = \(x) {
    n2 <- n2_extrapol(design, x[1L,])
    se2 <- sigma * sqrt((1L + two_armed) /  n2)
    smean1_prime <- z_to_smean(x[1L,], n1, sigma, two_armed)
    f1_kv(x[1L,,drop=FALSE], n1, mu, sigma, two_armed) *
      pnorm( ((smean * sqrtn * sqrt(n1 + n2) + mu * (n1 + n2 - sqrtn*sqrt(n1 + n2)) - n1 *  smean1_prime ) / n2)  /se2,
             mean = mu/se2,lower.tail = FALSE)
  },
  lowerLimit = c(design@c1f),
  upperLimit = c(design@c1e),
  tol = tol,
  maxEval = maxEval,
  absError = absError,
  vectorInterface = TRUE)$integral
  return(p_part_first_stage + p_part_second_stage)
}
p1_lr <- function(design, smean1, n1, mu, sigma, two_armed, ...) p_lr(design, smean1, n1, mu, sigma, two_armed, ...)
p2_lr <- function(design, smean1, smean2, n1, n2, mu, sigma, two_armed, ...) p_lr(design, smeans_to_smean(smean1, smean2, n1, n2), (n1 + n2), mu, sigma, two_armed, ...)

p_st <- function(design, smean, n, mu, sigma, two_armed, tol = getOption("adestr_tol_inner", default = .adestr_options[["adestr_tol_inner"]]), maxEval = getOption("adestr_maxEval_inner", default = .adestr_options[["adestr_maxEval_inner"]]), absError = getOption("adestr_absError_inner", default = .adestr_options[["adestr_absError_inner"]]), ...) {
  design <- TwoStageDesignWithCache(design)
  n1 <- n1(design, round=FALSE)
  se1 <- sigma_to_se(sigma, n1, two_armed)
  muse1 <- mu/se1
  smean1_sol <- (smean * n + mu*n1 - mu*n)/n1
  z1 <- smean_to_z(smean1_sol, n1, sigma, two_armed)
  if (z1 < design@c1f) {
    p_part_first_stage <- pnorm(design@c1f, mean = muse1) - pnorm(z1, mean = muse1) + pnorm(design@c1e, mean = muse1, lower.tail = FALSE)
  } else if (z1 > design@c1e) {
    p_part_first_stage <- pnorm(z1, mean = muse1, lower.tail = FALSE)
  } else {
    p_part_first_stage <- pnorm(design@c1e, mean = muse1, lower.tail = FALSE)
  }
  p_part_second_stage <- .hcubature(f = \(x) {
    n2 <- n2_extrapol(design, x[1L,])
    se2 <- sigma * sqrt((1L + two_armed) /  n2)
    smean1_prime <- z_to_smean(x[1L,], n1, sigma, two_armed)
    f1_kv(x[1L,,drop=FALSE], n1, mu, sigma, two_armed) *
      pnorm( ((smean * n + mu * (n1 + n2 - n) - n1*smean1_prime)/n2)  /se2,
             mean = mu/se2,lower.tail = FALSE)
  },
  lowerLimit = c(design@c1f),
  upperLimit = c(design@c1e),
  tol = tol,
  maxEval = maxEval,
  absError = absError,
  vectorInterface = TRUE)$integral
  return(p_part_first_stage + p_part_second_stage)
}
p1_st <- function(design, smean1, n1, mu, sigma, two_armed, ...) p_st(design, smean1, n1, mu, sigma, two_armed, ...)
p2_st <- function(design, smean1, smean2, n1, n2, mu, sigma, two_armed, ...) p_st(design, smeans_to_smean(smean1, smean2, n1, n2), (n1 + n2), mu, sigma, two_armed, ...)

p1_sw <- function(smean1, n1, mu, sigma, two_armed, ...) {
  se1 <- sigma_to_se(sigma, n1, two_armed)
  return(pnorm(smean_to_z(smean1, n1, sigma, two_armed), mean = mu/se1, lower.tail=FALSE))
}
p2_sw <- function(design, smean1, smean2, n1, n2, mu, sigma, two_armed, tol = getOption("adestr_tol_inner", default = .adestr_options[["adestr_tol_inner"]]), maxEval = getOption("adestr_maxEval_inner", default = .adestr_options[["adestr_maxEval_inner"]]), absError = getOption("adestr_absError_inner", default = .adestr_options[["adestr_absError_inner"]]), ...) {
  design <- TwoStageDesignWithCache(design)
  se1 <- sigma_to_se(sigma, n1, two_armed)
  CZ1Z2_shifted <- combf_kv(design, smean1, smean2, n1, n2, mu, sigma, two_armed)
  return(.hcubature(f = \(x) {
    f1_kv(x[1L,,drop=FALSE], n1, mu, sigma, two_armed) *
      pnorm( CZ1Z2_shifted + c2_extrapol(design, x[1L, ] - mu / se1),
             lower.tail = FALSE)
  },
  lowerLimit = c(design@c1f),
  upperLimit = c(design@c1e),
  tol = tol,
  maxEval = maxEval,
  absError = absError,
  vectorInterface = TRUE)$integral + pnorm(design@c1e, mean = mu/se1, lower.tail = FALSE))
}
p_sw <- function(design, smean1, smean2, n1, n2, mu, sigma, two_armed, tol = getOption("adestr_tol_inner", default = .adestr_options[["adestr_tol_inner"]]), maxEval = getOption("adestr_maxEval_inner", default = .adestr_options[["adestr_maxEval_inner"]]), absError = getOption("adestr_absError_inner", default = .adestr_options[["adestr_absError_inner"]]), ...) {
  z1 <- smean_to_z(smean1, n1, sigma, two_armed)
  if (z1 < design@c1f || z1 > design@c1f) {
    p1_sw(smean1, n1, mu, sigma, two_armed, ...)
  } else {
    p2_sw(design, smean1, smean2, n1, n2, mu, sigma, two_armed, tol, ...)
  }
}
find_root_p1_sw <- function(smean1, n1, sigma, two_armed, p_boundary, ...) {
  se1 <- sigma_to_se(sigma, n1, two_armed)
  return(qnorm(p_boundary, mean = smean1, sd = se1))
}
find_root_p2_sw <- function(design, smean1, smean2, n1, n2, sigma, two_armed, p_boundary, tol = getOption("adestr_tol_roots", default = .adestr_options[["adestr_tol_roots"]]), maxiter = getOption("adestr_maxiter_roots", default = .adestr_options[["adestr_maxiter_roots"]]), ...) {
  sec <- sigma_to_se(sigma, n1 + n2, two_armed)
  initial_guess <- (smean1*n1 + smean2*n2)/ (n1 + n2) + qnorm(p_boundary, sd = sec)
  uniroot(f = \(x) p2_sw(design = design, smean1 = smean1, smean2 = smean2,
                         n1 = n1, n2 = n2, mu = x, sigma = sigma, two_armed = two_armed,
                         tol = getOption("adestr_tol_inner", default = .adestr_options[["adestr_tol_inner"]]), maxEval = getOption("adestr_maxEval_inner", default = .adestr_options[["adestr_maxEval_inner"]]), absError = getOption("adestr_absError_inner", default = .adestr_options[["adestr_absError_inner"]])) - p_boundary,
          interval = c(initial_guess - sec,
                       initial_guess + sec),
          tol = tol,
          maxiter = maxiter,
          extendInt = "yes")$root
}

find_root_p1_ml <- function(design, smean1, n1, sigma, two_armed, p_boundary, tol = getOption("adestr_tol_roots", default = .adestr_options[["adestr_tol_roots"]]), maxiter = getOption("adestr_maxiter_roots", default = .adestr_options[["adestr_maxiter_roots"]]), ...) {
  se1 <- sigma_to_se(sigma, n1, two_armed)
  initial_guess <- smean1 + qnorm(p_boundary, sd = se1)
  uniroot(f = \(x) p1_ml(design = design, smean1 = smean1,
                         n1 = n1, mu = x, sigma = sigma, two_armed = two_armed,
                         tol = getOption("adestr_tol_inner", default = .adestr_options[["adestr_tol_inner"]]), maxEval = getOption("adestr_maxEval_inner", default = .adestr_options[["adestr_maxEval_inner"]]), absError = getOption("adestr_absError_inner", default = .adestr_options[["adestr_absError_inner"]])) - p_boundary,
          interval = c(initial_guess - se1,
                       initial_guess + se1),
          tol = tol,
          maxiter = maxiter,
          extendInt = "yes")$root
}
find_root_p2_ml <- function(design, smean1, smean2, n1, n2, sigma, two_armed, p_boundary, tol = getOption("adestr_tol_roots", default = .adestr_options[["adestr_tol_roots"]]), maxiter = getOption("adestr_maxiter_roots", default = .adestr_options[["adestr_maxiter_roots"]]), ...) {
  se <- sigma_to_se(sigma, n1 + n2, two_armed)
  initial_guess <- smeans_to_smean(smean1, smean2, n1, n2) + qnorm(p_boundary, sd = se)
  uniroot(f = \(x) p2_ml(design = design, smean1 = smean1, smean2 = smean2,
                         n1 = n1, n2 = n2, mu = x, sigma = sigma, two_armed = two_armed,
                         tol = getOption("adestr_tol_inner", default = .adestr_options[["adestr_tol_inner"]]), maxEval = getOption("adestr_maxEval_inner", default = .adestr_options[["adestr_maxEval_inner"]]), absError = getOption("adestr_absError_inner", default = .adestr_options[["adestr_absError_inner"]])) - p_boundary,
          interval = c(initial_guess - se,
                       initial_guess + se),
          tol = tol,
          maxiter = maxiter,
          extendInt = "yes")$root
}
find_root_p1_lr <- function(design, smean1, n1, sigma, two_armed, p_boundary, tol = getOption("adestr_tol_roots", default = .adestr_options[["adestr_tol_roots"]]), maxiter = getOption("adestr_maxiter_roots", default = .adestr_options[["adestr_maxiter_roots"]]), ...) {
  se1 <- sigma_to_se(sigma, n1, two_armed)
  initial_guess <- smean1 + qnorm(p_boundary, sd = se1)
  uniroot(f = \(x) p1_lr(design = design, smean1 = smean1,
                         n1 = n1, mu = x, sigma = sigma, two_armed = two_armed,
                         tol = getOption("adestr_tol_inner", default = .adestr_options[["adestr_tol_inner"]]), maxEval = getOption("adestr_maxEval_inner", default = .adestr_options[["adestr_maxEval_inner"]]), absError = getOption("adestr_absError_inner", default = .adestr_options[["adestr_absError_inner"]])) - p_boundary,
          interval = c(initial_guess - se1,
                       initial_guess + se1),
          tol = tol,
          maxiter = maxiter,
          extendInt = "yes")$root
}
find_root_p2_lr <- function(design, smean1, smean2, n1, n2, sigma, two_armed, p_boundary, tol = getOption("adestr_tol_roots", default = .adestr_options[["adestr_tol_roots"]]), maxiter = getOption("adestr_maxiter_roots", default = .adestr_options[["adestr_maxiter_roots"]]), ...) {
  se <- sigma_to_se(sigma, n1 + n2, two_armed)
  initial_guess <- smeans_to_smean(smean1, smean2, n1, n2) + qnorm(p_boundary, sd = se)
  uniroot(f = \(x) p2_lr(design = design, smean1 = smean1, smean2 = smean2,
                         n1 = n1, n2 = n2, mu = x, sigma = sigma, two_armed = two_armed,
                         tol = getOption("adestr_tol_inner", default = .adestr_options[["adestr_tol_inner"]]), maxEval = getOption("adestr_maxEval_inner", default = .adestr_options[["adestr_maxEval_inner"]]), absError = getOption("adestr_absError_inner", default = .adestr_options[["adestr_absError_inner"]])) - p_boundary,
          interval = c(initial_guess - se,
                       initial_guess + se),
          tol = tol,
          maxiter = maxiter,
          extendInt = "yes")$root
}
find_root_p1_st <- function(design, smean1, n1, sigma, two_armed, p_boundary, tol = getOption("adestr_tol_roots", default = .adestr_options[["adestr_tol_roots"]]), maxiter = getOption("adestr_maxiter_roots", default = .adestr_options[["adestr_maxiter_roots"]]), ...) {
  se1 <- sigma_to_se(sigma, n1, two_armed)
  initial_guess <- smean1 + qnorm(p_boundary, sd = se1)
  uniroot(f = \(x) p1_st(design = design, smean1 = smean1,
                         n1 = n1, mu = x, sigma = sigma, two_armed = two_armed,
                         tol = getOption("adestr_tol_inner", default = .adestr_options[["adestr_tol_inner"]]), maxEval = getOption("adestr_maxEval_inner", default = .adestr_options[["adestr_maxEval_inner"]]), absError = getOption("adestr_absError_inner", default = .adestr_options[["adestr_absError_inner"]])) - p_boundary,
          interval = c(initial_guess + 2.5 * se1,
                       initial_guess - 2.5 * se1),
          tol = tol,
          maxiter = maxiter,
          extendInt = "yes")$root
}
find_root_p2_st <- function(design, smean1, smean2, n1, n2, sigma, two_armed, p_boundary, tol = getOption("adestr_tol_roots", default = .adestr_options[["adestr_tol_roots"]]), maxiter = getOption("adestr_maxiter_roots", default = .adestr_options[["adestr_maxiter_roots"]]), ...) {
  se <- sigma_to_se(sigma, n1 + n2, two_armed)
  initial_guess <- smeans_to_smean(smean1, smean2, n1, n2) + qnorm(p_boundary, sd = se)
  uniroot(f = \(x) p2_st(design = design, smean1 = smean1, smean2 = smean2,
                         n1 = n1, n2 = n2, mu = x, sigma = sigma, two_armed = two_armed,
                         tol = getOption("adestr_tol_inner", default = .adestr_options[["adestr_tol_inner"]]), maxEval = getOption("adestr_maxEval_inner", default = .adestr_options[["adestr_maxEval_inner"]]), absError = getOption("adestr_absError_inner", default = .adestr_options[["adestr_absError_inner"]])) - p_boundary,
          interval = c(initial_guess - se,
                       initial_guess + se),
          tol = tol,
          maxiter = maxiter,
          extendInt = "yes")$root
}
p_np <- function(design, smean, n, mu, mu0, mu1, sigma, two_armed, tol = getOption("adestr_tol_inner", default = .adestr_options[["adestr_tol_inner"]]), maxEval = getOption("adestr_maxEval_inner", default = .adestr_options[["adestr_maxEval_inner"]]), absError = getOption("adestr_absError_inner", default = .adestr_options[["adestr_absError_inner"]]), ...) {
  design <- TwoStageDesignWithCache(design)
  n1 <- n1(design, round=FALSE)
  se1 <- sigma_to_se(sigma, n1, two_armed)
  muse1 <- mu/se1
  smean1_sol <- (2L*smean*n + mu0*(n1 - n) + mu1*(n1 - n))/(2L*n1)
  z1 <- smean_to_z(smean1_sol, n1, sigma, two_armed)
  if (z1 < design@c1f) {
    p_part_first_stage <- pnorm(design@c1f, mean = muse1) - pnorm(z1, mean = muse1) + pnorm(design@c1e, mean = muse1, lower.tail = FALSE)
  } else if (z1 > design@c1e) {
    p_part_first_stage <- pnorm(z1, mean = muse1, lower.tail = FALSE)
  } else {
    p_part_first_stage <- pnorm(design@c1e, mean = muse1, lower.tail = FALSE)
  }
  const <- n*(2*smean - mu0 - mu1)

  rr <- .rho_np2(x1 = .2, x2 = .2, n1 = n1, n2 = n-n1, mu0 = mu0, mu1 = mu1)
  p_part_second_stage <- .hcubature(f = \(x) {
    n2 <- n2_extrapol(design, x[1L,])
    se2 <- sigma * sqrt((1L + two_armed) /  n2)
    smean1_prime <- z_to_smean(x[1L,], n1, sigma, two_armed)
    nn <- (n1 + n2)
    f1_kv(x[1L,,drop=FALSE], n1, mu, sigma, two_armed) *
      pnorm( (const + mu0*nn + mu1*nn - 2*n1*smean1_prime)  / (2 * n2 * se2),
             mean = mu/se2,lower.tail = FALSE)
  },
  lowerLimit = c(design@c1f),
  upperLimit = c(design@c1e),
  tol = tol,
  maxEval = maxEval,
  absError = absError,
  vectorInterface = TRUE)$integral
  return(p_part_first_stage + p_part_second_stage)
}
p1_np <- function(design, smean1, n1, mu, mu0, mu1, sigma, two_armed, ...) p_np(design, smean1, n1, mu, mu0, mu1, sigma, two_armed, ...)
p2_np <- function(design, smean1, smean2, n1, n2, mu, mu0, mu1, sigma, two_armed, ...) p_np(design, smeans_to_smean(smean1, smean2, n1, n2), n1+n2, mu, mu0, mu1, sigma, two_armed, ...)

find_root_p1_np <- function(design, smean1, n1, mu0, mu1, sigma, two_armed, p_boundary, tol = getOption("adestr_tol_roots", default = .adestr_options[["adestr_tol_roots"]]), maxiter = getOption("adestr_maxiter_roots", default = .adestr_options[["adestr_maxiter_roots"]]), ...) {
  se1 <- sigma_to_se(sigma, n1, two_armed)
  initial_guess <- smean1 + qnorm(p_boundary, sd = se1)
  uniroot(f = \(x) p1_np(design = design, smean1 = smean1,
                         n1 = n1, mu = x, mu0 = mu0 +x, mu1 = mu1 +x, sigma = sigma, two_armed = two_armed,
                         tol = getOption("adestr_tol_inner", default = .adestr_options[["adestr_tol_inner"]]), maxEval = getOption("adestr_maxEval_inner", default = .adestr_options[["adestr_maxEval_inner"]]), absError = getOption("adestr_absError_inner", default = .adestr_options[["adestr_absError_inner"]])) - p_boundary,
          interval = c(initial_guess - se1,
                       initial_guess + se1),
          tol = tol,
          maxiter = maxiter,
          extendInt = "yes")$root
}
find_root_p2_np <- function(design, smean1, smean2, n1, n2, mu0, mu1, sigma, two_armed, p_boundary, tol = getOption("adestr_tol_roots", default = .adestr_options[["adestr_tol_roots"]]), maxiter = getOption("adestr_maxiter_roots", default = .adestr_options[["adestr_maxiter_roots"]]), ...) {
  se <- sigma_to_se(sigma, n1 + n2, two_armed)
  initial_guess <- smeans_to_smean(smean1, smean2, n1, n2) + qnorm(p_boundary, sd = se)
  uniroot(f = \(x) p2_np(design = design, smean1 = smean1, smean2 = smean2,
                         n1 = n1, n2 = n2, mu = x, mu0 = mu0 + x, mu1 = mu1 + x, sigma = sigma, two_armed = two_armed,
                         tol = getOption("adestr_tol_inner", default = .adestr_options[["adestr_tol_inner"]]), maxEval = getOption("adestr_maxEval_inner", default = .adestr_options[["adestr_maxEval_inner"]]), absError = getOption("adestr_absError_inner", default = .adestr_options[["adestr_absError_inner"]])) - p_boundary,
          interval = c(initial_guess - se,
                       initial_guess + se),
          tol = tol,
          maxiter = maxiter,
          extendInt = "yes")$root
}

setClass("MLEOrderingPValue", contains = "VirtualPValue")
#' @rdname PValue-class
#' @export
MLEOrderingPValue <- function() new("MLEOrderingPValue", label = "MLE ordering p-value")
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("MLEOrderingPValue", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            g1 <- Vectorize(\(design, smean1, n1, sigma, two_armed, ...) p1_ml(design, smean1, n1, mu = 0, sigma, two_armed, ...), c("smean1"))
            g2 <- Vectorize(\(design, smean1, smean2, n1, n2, sigma, two_armed, ...) p2_ml(design, smean1, smean2, n1, n2, mu = 0, sigma, two_armed, ...), c("smean1", "smean2", "n2"))
            list(g1 = g1,
                 g2 = g2)
          })
setClass("LikelihoodRatioOrderingPValue", contains = "VirtualPValue")
#' @rdname PValue-class
#' @export
LikelihoodRatioOrderingPValue <- function() new("LikelihoodRatioOrderingPValue", label = "LR test ordering p-value")
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("LikelihoodRatioOrderingPValue", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            g1 <- Vectorize(\(design, smean1, n1, sigma, two_armed, ...) p1_lr(design, smean1, n1, mu = 0, sigma, two_armed, ...), c("smean1"))
            g2 <- Vectorize(\(design, smean1, smean2, n1, n2, sigma, two_armed, ...) p2_lr(design, smean1, smean2, n1, n2, mu = 0, sigma, two_armed, ...), c("smean1", "smean2", "n2"))
            list(g1 = g1,
                 g2 = g2)
          })
setClass("ScoreTestOrderingPValue", contains = "VirtualPValue")
#' @rdname PValue-class
#' @export
ScoreTestOrderingPValue <- function() new("ScoreTestOrderingPValue", label = "Score test ordering p-value")
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("ScoreTestOrderingPValue", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            g1 <- Vectorize(\(smean1, n1, sigma, two_armed, ...) p1_st(smean1, n1, mu = 0, sigma, two_armed, ...), c("smean1"))
            g2 <- Vectorize(\(design, smean1, smean2, n1, n2, sigma, two_armed, ...) p2_st(design, smean1, smean2, n1, n2, mu = 0, sigma, two_armed, ...), c("smean1", "smean2", "n2"))
            list(g1 = g1,
                 g2 = g2)
          })
setClass("StagewiseCombinationFunctionOrderingPValue", contains = "VirtualPValue")
#' @rdname PValue-class
#' @export
StagewiseCombinationFunctionOrderingPValue <- function() new("StagewiseCombinationFunctionOrderingPValue", label = "SWCF ordering p-value")
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("StagewiseCombinationFunctionOrderingPValue", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            g1 <- Vectorize(\(smean1, n1, sigma, two_armed, ...) p1_sw(smean1, n1, mu = 0, sigma, two_armed, ...), c("smean1"))
            g2 <- Vectorize(\(design, smean1, smean2, n1, n2, sigma, two_armed, ...) p2_sw(design, smean1, smean2, n1, n2, mu = 0, sigma, two_armed, ...), c("smean1", "smean2", "n2"))
            list(g1 = g1,
                 g2 = g2)
          })
setClass("NeymanPearsonOrderingPValue", contains = "VirtualPValue", slots = c(mu0 = "numeric", mu1 = "numeric"))
#' @param mu0 expected value of the normal distribution under the null hypothesis.
#' @param mu1 expected value of the normal distribution under the null hypothesis.
#' @rdname PValue-class
#' @export
NeymanPearsonOrderingPValue <- function(mu0 = 0, mu1 = 0.4) new("NeymanPearsonOrderingPValue", mu0 =  mu0, mu1 = mu1, label = paste0("Neyman-Pearson test ordering (mu0=",format(mu0),", mu1=",format(mu1),")") )
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("NeymanPearsonOrderingPValue", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            g1 <- Vectorize(\(design, smean1, n1, sigma, two_armed, ...) p1_np(design, smean1, n1, mu = 0, mu0=estimator@mu0, mu1=estimator@mu1, sigma, two_armed, ...), c("smean1"))
            g2 <- Vectorize(\(design, smean1, smean2, n1, n2, sigma, two_armed, ...) p2_np(design, smean1, smean2, n1, n2, mu = 0, mu0=estimator@mu0, mu1=estimator@mu1, sigma, two_armed, ...), c("smean1", "smean2", "n2"))
            list(g1 = g1,
                 g2 = g2)
          })

setClass("StagewiseCombinationFunctionOrderingCI", contains = "VirtualIntervalEstimator")
#' @rdname IntervalEstimator-class
#' @export
StagewiseCombinationFunctionOrderingCI <- function(two_sided = TRUE) new("StagewiseCombinationFunctionOrderingCI", two_sided = two_sided, label = "SWCF ordering CI")
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("StagewiseCombinationFunctionOrderingCI", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            alpha <- adoptr_alpha_shifted_design_kv(design, 0, 0, 0)
            l1 <- Vectorize(\(smean1, n1, sigma, two_armed, ...) find_root_p1_sw(smean1, n1, sigma, two_armed, p_boundary = alpha, ...), c("smean1"))
            l2 <- Vectorize(\(design, smean1, smean2, n1, n2, sigma, two_armed, ...) find_root_p2_sw(design, smean1, smean2, n1, n2, sigma, two_armed, p_boundary = alpha, ...), c("smean1", "smean2", "n2"))
            if (estimator@two_sided){
              u1 <- Vectorize(\(smean1, n1, sigma, two_armed, ...) find_root_p1_sw(smean1, n1, sigma, two_armed, p_boundary = 1-alpha, ...), c("smean1"))
              u2 <- Vectorize(\(design, smean1, smean2, n1, n2, sigma, two_armed, ...) find_root_p2_sw(design, smean1, smean2, n1, n2, sigma, two_armed, p_boundary = 1-alpha, ...), c("smean1", "smean2", "n2"))
            } else {
              u1 <- \(...) Inf
              u2 <- \(...) Inf
            }
            return(list(
              l1 = l1,
              u1 = u1,
              l2 = l2,
              u2 = u2
            ))
          })

setClass("MLEOrderingCI", contains = "VirtualIntervalEstimator")
#' @rdname IntervalEstimator-class
#' @export
MLEOrderingCI <- function(two_sided = TRUE) new("MLEOrderingCI", two_sided = two_sided, label = "MLE ordering CI")
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("MLEOrderingCI", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            alpha <- adoptr_alpha_shifted_design_kv(design, 0, 0, 0)
            l1 <- Vectorize(\(design, smean1, n1, sigma, two_armed, ...) find_root_p1_ml(design, smean1, n1, sigma, two_armed, p_boundary = alpha, ...), c("smean1"))
            l2 <- Vectorize(\(design, smean1, smean2, n1, n2, sigma, two_armed, ...) find_root_p2_ml(design, smean1, smean2, n1, n2, sigma, two_armed, p_boundary = alpha, ...), c("smean1", "smean2", "n2"))
            if (estimator@two_sided){
              u1 <- Vectorize(\(design, smean1, n1, sigma, two_armed, ...) find_root_p1_ml(design, smean1, n1, sigma, two_armed, p_boundary = 1-alpha, ...), c("smean1"))
              u2 <- Vectorize(\(design, smean1, smean2, n1, n2, sigma, two_armed, ...) find_root_p2_ml(design, smean1, smean2, n1, n2, sigma, two_armed, p_boundary = 1-alpha, ...), c("smean1", "smean2", "n2"))
            } else {
              u1 <- \(design, smean1, sigma, two_armed, ...) Inf
              u2 <- \(design, smean1, smean2, sigma, two_armed, ...) Inf
            }
            return(list(
              l1 = l1,
              u1 = u1,
              l2 = l2,
              u2 = u2
            ))
          })
setClass("LikelihoodRatioOrderingCI", contains = "VirtualIntervalEstimator")
#' @rdname IntervalEstimator-class
#' @export
LikelihoodRatioOrderingCI <- function(two_sided = TRUE) new("LikelihoodRatioOrderingCI", two_sided = two_sided, label = "LR test ordering CI")
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("LikelihoodRatioOrderingCI", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            alpha <- adoptr_alpha_shifted_design_kv(design, 0, 0, 0)
            l1 <- Vectorize(\(design, smean1, n1, sigma, two_armed, ...) find_root_p1_lr(design, smean1, n1, sigma, two_armed, p_boundary = alpha, ...), c("smean1"))
            l2 <- Vectorize(\(design, smean1, smean2, n1, n2, sigma, two_armed, ...) find_root_p2_lr(design, smean1, smean2, n1, n2, sigma, two_armed, p_boundary = alpha, ...), c("smean1", "smean2", "n2"))
            if (estimator@two_sided){
              u1 <- Vectorize(\(design, smean1, n1, sigma, two_armed, ...) find_root_p1_lr(design, smean1, n1, sigma, two_armed, p_boundary = 1-alpha, ...), c("smean1"))
              u2 <- Vectorize(\(design, smean1, smean2, n1, n2, sigma, two_armed, ...) find_root_p2_lr(design, smean1, smean2, n1, n2, sigma, two_armed, p_boundary = 1-alpha, ...), c("smean1", "smean2", "n2"))
            } else {
              u1 <- \(design, smean1, sigma, two_armed, ...) Inf
              u2 <- \(design, smean1, smean2, sigma, two_armed, ...) Inf
            }
            return(list(
              l1 = l1,
              u1 = u1,
              l2 = l2,
              u2 = u2
            ))
          })
setClass("ScoreTestOrderingCI", contains = "VirtualIntervalEstimator")
#' @rdname IntervalEstimator-class
#' @export
ScoreTestOrderingCI <- function(two_sided = TRUE) new("ScoreTestOrderingCI", two_sided = two_sided, label = "Score test ordering CI")
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("ScoreTestOrderingCI", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            alpha <- adoptr_alpha_shifted_design_kv(design, 0, 0, 0)
            l1 <- Vectorize(\(design, smean1, n1, sigma, two_armed, ...) find_root_p1_st(design, smean1, n1, sigma, two_armed, p_boundary = alpha, ...), c("smean1"))
            l2 <- Vectorize(\(design, smean1, smean2, n1, n2, sigma, two_armed, ...) find_root_p2_st(design, smean1, smean2, n1, n2, sigma, two_armed, p_boundary = alpha, ...), c("smean1", "smean2", "n2"))
            if (estimator@two_sided){
              u1 <- Vectorize(\(design, smean1, n1, sigma, two_armed, ...) find_root_p1_st(design, smean1, n1, sigma, two_armed, p_boundary = 1-alpha, ...), c("smean1"))
              u2 <- Vectorize(\(design, smean1, smean2, n1, n2, sigma, two_armed, ...) find_root_p2_st(design, smean1, smean2, n1, n2, sigma, two_armed, p_boundary = 1-alpha, ...), c("smean1", "smean2", "n2"))
            } else {
              u1 <- \(design, smean1, sigma, two_armed, ...) Inf
              u2 <- \(design, smean1, smean2, sigma, two_armed, ...) Inf
            }
            return(list(
              l1 = l1,
              u1 = u1,
              l2 = l2,
              u2 = u2
            ))
          })
setClass("NeymanPearsonOrderingCI", contains = "VirtualIntervalEstimator", slots = c(mu0 = "numeric", mu1 = "numeric"))
#' @inheritParams NeymanPearsonOrderingPValue
#' @rdname IntervalEstimator-class
#' @export
NeymanPearsonOrderingCI <- function(two_sided = TRUE, mu0 = 0, mu1 = 0.4) new("NeymanPearsonOrderingCI", mu0 = mu0, mu1 = mu1, two_sided = two_sided, label = paste0("Neyman-Pearson test ordering (mu0=",format(mu0),", mu1=",format(mu1),")"))
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("NeymanPearsonOrderingCI", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            alpha <- adoptr_alpha_shifted_design_kv(design, 0, 0, 0)
            l1 <- Vectorize(\(design, smean1, n1, sigma, two_armed, ...) find_root_p1_np(design, smean1, n1, mu0=estimator@mu0, mu1=estimator@mu1, sigma, two_armed, p_boundary = alpha, ...), c("smean1"))
            l2 <- Vectorize(\(design, smean1, smean2, n1, n2, sigma, two_armed, ...) find_root_p2_np(design, smean1, smean2, n1, n2, mu0=estimator@mu0, mu1=estimator@mu1, sigma, two_armed, p_boundary = alpha, ...), c("smean1", "smean2", "n2"))
            if (estimator@two_sided){
              u1 <- Vectorize(\(design, smean1, n1, sigma, two_armed, ...) find_root_p1_np(design, smean1, n1, mu0=estimator@mu0, mu1=estimator@mu1, sigma, two_armed, p_boundary = 1-alpha, ...), c("smean1"))
              u2 <- Vectorize(\(design, smean1, smean2, n1, n2, sigma, two_armed, ...) find_root_p2_np(design, smean1, smean2, n1, n2, mu0=estimator@mu0, mu1=estimator@mu1, sigma, two_armed, p_boundary = 1-alpha, ...), c("smean1", "smean2", "n2"))
            } else {
              u1 <- \(design, smean1, sigma, two_armed, ...) Inf
              u2 <- \(design, smean1, smean2, sigma, two_armed, ...) Inf
            }
            return(list(
              l1 = l1,
              u1 = u1,
              l2 = l2,
              u2 = u2
            ))
          })

setClass("NaiveCI", contains = "VirtualIntervalEstimator")
#' @rdname IntervalEstimator-class
#' @export
NaiveCI <- function(two_sided = TRUE) new("NaiveCI", two_sided = two_sided, label = "Naive CI")
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("NaiveCI", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            alpha <- adoptr_alpha_shifted_design_kv(design, 0, 0, 0)
            l1 <- \(smean1, n1, sigma, two_armed, ...) {
              qnorm(alpha, mean = smean1, sd = sigma * sqrt((1L + two_armed) / n1), lower.tail = TRUE)
            }
            l2 <- \(smean1, smean2, n1, n2, sigma, two_armed, ...) {
              qnorm(alpha, mean = (smean1 * n1 + smean2 * n2) / (n1 + n2), sd = sigma * sqrt((1L + two_armed) / (n1 + n2)), lower.tail = TRUE)
            }
            if (estimator@two_sided){
              u1 <- \(smean1, n1, sigma, two_armed, ...) {
                qnorm(alpha, mean = smean1, sd = sigma * sqrt((1L + two_armed) / n1), lower.tail = FALSE)
              }
              u2 <- \(smean1, smean2, n1, n2, sigma, two_armed, ...) {
                qnorm(alpha, mean = (smean1 * n1 + smean2 * n2) / (n1 + n2), sd = sigma * sqrt((1L + two_armed) / (n1 + n2)), lower.tail = FALSE)
              }
            } else {
              u1 <- \(smean1, n1, sigma, two_armed, ...) {
                Inf
              }
              u2 <- \(smean1, smean2, n1, n2, sigma, two_armed, ...) {
                Inf
              }
            }
            return(list(
              l1 = l1,
              u1 = u1,
              l2 = l2,
              u2 = u2
            ))
          })

setClass("MidpointStagewiseCombinationFunctionOrderingCI", contains = "VirtualPointEstimator")
#' @rdname PointEstimator-class
#' @export
MidpointStagewiseCombinationFunctionOrderingCI <- function() new("MidpointStagewiseCombinationFunctionOrderingCI", label = "Midpoint of SWCF ordering CI")
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("MidpointStagewiseCombinationFunctionOrderingCI", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            cis <- get_stagewise_estimators(estimator = StagewiseCombinationFunctionOrderingCI(two_sided = TRUE),
                                            data_distribution =  data_distribution,
                                            use_full_twoarm_sampling_distribution = use_full_twoarm_sampling_distribution,
                                            design = design, sigma = sigma, exact = exact)
            l1 <- cis$l1
            l2 <- cis$l2
            u1 <- cis$u1
            u2 <- cis$u2
            g1 <- \(...) (l1(...) + u1(...))/2
            g2 <- \(...) (l2(...) + u2(...))/2
            return(list(
              g1 = g1,
              g2 = g2
            ))
          })
setClass("MidpointMLEOrderingCI", contains = "VirtualPointEstimator")
#' @rdname PointEstimator-class
#' @export
MidpointMLEOrderingCI <- function() new("MidpointMLEOrderingCI", label = "Midpoint of MLE ordering CI")
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("MidpointMLEOrderingCI", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            cis <- get_stagewise_estimators(estimator = MLEOrderingCI(two_sided = TRUE),
                                            data_distribution =  data_distribution,
                                            use_full_twoarm_sampling_distribution = use_full_twoarm_sampling_distribution,
                                            design = design, sigma = sigma, exact = exact)
            l1 <- cis$l1
            l2 <- cis$l2
            u1 <- cis$u1
            u2 <- cis$u2
            g1 <- \(...) (l1(...) + u1(...))/2
            g2 <- \(...) (l2(...) + u2(...))/2
            return(list(
              g1 = g1,
              g2 = g2
            ))
          })
setClass("MidpointLikelihoodRatioOrderingCI", contains = "VirtualPointEstimator")
#' @rdname PointEstimator-class
#' @export
MidpointLikelihoodRatioOrderingCI <- function() new("MidpointLikelihoodRatioOrderingCI", label = "Midpoint of LR test ordering CI")
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("MidpointLikelihoodRatioOrderingCI", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            cis <- get_stagewise_estimators(estimator = LikelihoodRatioOrderingCI(two_sided = TRUE),
                                            data_distribution =  data_distribution,
                                            use_full_twoarm_sampling_distribution = use_full_twoarm_sampling_distribution,
                                            design = design, sigma = sigma, exact = exact)
            l1 <- cis$l1
            l2 <- cis$l2
            u1 <- cis$u1
            u2 <- cis$u2
            g1 <- \(...) (l1(...) + u1(...))/2
            g2 <- \(...) (l2(...) + u2(...))/2
            return(list(
              g1 = g1,
              g2 = g2
            ))
          })
setClass("MidpointScoreTestOrderingCI", contains = "VirtualPointEstimator")
#' @rdname PointEstimator-class
#' @export
MidpointScoreTestOrderingCI <- function() new("MidpointScoreTestOrderingCI", label = "Midpoint of score test ordering CI")
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("MidpointScoreTestOrderingCI", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            cis <- get_stagewise_estimators(estimator = ScoreTestOrderingCI(two_sided = TRUE),
                                            data_distribution =  data_distribution,
                                            use_full_twoarm_sampling_distribution = use_full_twoarm_sampling_distribution,
                                            design = design, sigma = sigma, exact = exact)
            l1 <- cis$l1
            l2 <- cis$l2
            u1 <- cis$u1
            u2 <- cis$u2
            g1 <- \(...) (l1(...) + u1(...))/2
            g2 <- \(...) (l2(...) + u2(...))/2
            return(list(
              g1 = g1,
              g2 = g2
            ))
          })
setClass("MidpointNeymanPearsonOrderingCI", contains = "VirtualPointEstimator")
#' @inheritParams NeymanPearsonOrderingPValue
#' @rdname PointEstimator-class
#' @export
MidpointNeymanPearsonOrderingCI <- function() new("MidpointNeymanPearsonOrderingCI", label = "Midpoint of Neyman-Pearson ordering CI")
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("MidpointNeymanPearsonOrderingCI", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            cis <- get_stagewise_estimators(estimator = NeymanPearsonOrderingCI(two_sided = TRUE),
                                            data_distribution =  data_distribution,
                                            use_full_twoarm_sampling_distribution = use_full_twoarm_sampling_distribution,
                                            design = design, sigma = sigma, exact = exact)
            l1 <- cis$l1
            l2 <- cis$l2
            u1 <- cis$u1
            u2 <- cis$u2
            g1 <- \(...) (l1(...) + u1(...))/2
            g2 <- \(...) (l2(...) + u2(...))/2
            return(list(
              g1 = g1,
              g2 = g2
            ))
          })

setClass("MedianUnbiasedStagewiseCombinationFunctionOrdering", contains = "VirtualPointEstimator")
#' @rdname PointEstimator-class
#' @export
MedianUnbiasedStagewiseCombinationFunctionOrdering <- function() new("MedianUnbiasedStagewiseCombinationFunctionOrdering", label = "Median unbiased (SWCF ordering)")
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("MedianUnbiasedStagewiseCombinationFunctionOrdering", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            list(g1 = Vectorize(\(smean1, n1, sigma, two_armed, ...) find_root_p1_sw(smean1, n1, sigma, two_armed, p_boundary = 0.5, ...), c("smean1")),
                 g2 = Vectorize(\(design, smean1, smean2, n1, n2, sigma, two_armed, ...) find_root_p2_sw(design, smean1, smean2, n1, n2, sigma, two_armed, p_boundary = 0.5, ...), c("smean1", "smean2", "n2")))
          })
setClass("MedianUnbiasedMLEOrdering", contains = "VirtualPointEstimator")
#' @rdname PointEstimator-class
#' @export
MedianUnbiasedMLEOrdering <- function() new("MedianUnbiasedMLEOrdering", label = "Median unbiased (MLE ordering)")
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("MedianUnbiasedMLEOrdering", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            list(g1 = Vectorize(\(design, smean1, n1, sigma, two_armed, ...) find_root_p1_ml(design, smean1, n1, sigma, two_armed, p_boundary = 0.5, ...), c("smean1")),
                 g2 = Vectorize(\(design, smean1, smean2, n1, n2, sigma, two_armed, ...) find_root_p2_ml(design, smean1, smean2, n1, n2, sigma, two_armed, p_boundary = 0.5, ...), c("smean1", "smean2", "n2")))
          })
setClass("MedianUnbiasedLikelihoodRatioOrdering", contains = "VirtualPointEstimator")
#' @rdname PointEstimator-class
#' @export
MedianUnbiasedLikelihoodRatioOrdering <- function() new("MedianUnbiasedLikelihoodRatioOrdering", label = "Median unbiased (LR test ordering)")
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("MedianUnbiasedLikelihoodRatioOrdering", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            list(g1 = Vectorize(\(design, smean1, n1, sigma, two_armed, ...) find_root_p1_lr(design, smean1, n1, sigma, two_armed, p_boundary = 0.5, ...), c("smean1")),
                 g2 = Vectorize(\(design, smean1, smean2, n1, n2, sigma, two_armed, ...) find_root_p2_lr(design, smean1, smean2, n1, n2, sigma, two_armed, p_boundary = 0.5, ...), c("smean1", "smean2", "n2")))
          })
setClass("MedianUnbiasedScoreTestOrdering", contains = "VirtualPointEstimator")
#' @rdname PointEstimator-class
#' @export
MedianUnbiasedScoreTestOrdering <- function() new("MedianUnbiasedScoreTestOrdering", label = "Median unbiased (Score test ordering)")
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("MedianUnbiasedScoreTestOrdering", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            list(g1 = Vectorize(\(design, smean1, n1, sigma, two_armed, ...) find_root_p1_st(design, smean1, n1, sigma, two_armed, p_boundary = 0.5, ...), c("smean1")),
                 g2 = Vectorize(\(design, smean1, smean2, n1, n2, sigma, two_armed, ...) find_root_p2_st(design, smean1, smean2, n1, n2, sigma, two_armed, p_boundary = 0.5, ...), c("smean1", "smean2", "n2")))
          })
setClass("MedianUnbiasedNeymanPearsonOrdering", contains = "VirtualPointEstimator", slots = c(mu0 = "numeric", mu1 = "numeric"))
#' @inheritParams NeymanPearsonOrderingPValue
#' @rdname PointEstimator-class
#' @export
MedianUnbiasedNeymanPearsonOrdering <- function(mu0 = 0, mu1 = 0.4) new("MedianUnbiasedNeymanPearsonOrdering", mu0 = mu0, mu1 = mu1, label = paste0("Median unbiased (Neyman-Pearson test ordering, mu0=",format(mu0),", mu1=",format(mu1),")"))
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("MedianUnbiasedNeymanPearsonOrdering", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            list(g1 = Vectorize(\(design, smean1, n1, sigma, two_armed, ...) find_root_p1_np(design, smean1, n1, mu0=estimator@mu0, mu1=estimator@mu1, sigma, two_armed, p_boundary = 0.5, ...), c("smean1")),
                 g2 = Vectorize(\(design, smean1, smean2, n1, n2, sigma, two_armed, ...) find_root_p2_np(design, smean1, smean2, n1, n2, mu0=estimator@mu0, mu1=estimator@mu1, sigma, two_armed, p_boundary = 0.5, ...), c("smean1", "smean2", "n2")))
          })

# For the p-value plots
implied_c2 <- function(design, z1, p2, sigma, two_armed, alpha, tol_root = getOption("adestr_tol_roots", default = .adestr_options[["adestr_tol_roots"]]), tol = getOption("adestr_tol_inner", default = .adestr_options[["adestr_tol_inner"]]), maxEval = getOption("adestr_maxEval_inner", default = .adestr_options[["adestr_maxEval_inner"]]), absError = getOption("adestr_absError_inner", default = .adestr_options[["adestr_absError_inner"]]), ...){
  design <- TwoStageDesignWithCache(design)
  n1 <- n1(design, round=FALSE)
  initialz <- qnorm(alpha/2, lower.tail = FALSE)
  z1_ord <- order(z1)
  z1 <- z1[z1_ord]
  res <- sapply(z1,
                \(z) {
                  n2 <- n2_extrapol(design, z)
                  n <- n1 + n2
                  smean1 <- z_to_smean(z, n1, sigma, two_armed)
                  root <- uniroot(\(y) {
                    p2(
                      design = design,
                      smean1 = smean1,
                      smean2 = z_to_smean(y, n2, sigma, two_armed),
                      n1 = n1, n2 = n2,
                      sigma = sigma, two_armed = two_armed,
                      tol = tol, maxEval = maxEval, absError = absError
                    ) - alpha
                  } , lower = initialz - 1, upper = initialz + 1, extendInt = "yes", tol = tol_root)$root
                  initialz <<- root
                  root
                })
  res <- res[order(z1_ord)]
  res
}

