#' adestr
#' @description Adaptive Design Estimation in R
#'
#' @details The flexibility of adaptive clinical trial designs can offer significant advantages over fixed designs.
#' If design parameters are chosen appropriately, adaptive designs can decrease the expected sample size and
#' save resources in cases with early evidence that the continuation of a trial may be futile. A known issue with adaptive
#' designs is that estimators appropriate in a single-stage fixed design, such as the maximum likelihood estimator, can be
#' biased due to the dependence structure in the data introduced by adaptivity. This problem affects point estimators as well
#' as confidence intervals and p-values. Regulatory agencies such as the EMA or FDA recognize this problem and urge researchers
#' that "the extent of bias should be evaluated, and estimates should be presented with appropriate cautions regarding
#' their interpretation". Over the years, various methods have been proposed to mitigate the bias introduced by adaptive designs.
#' However, estimators often need to fulfill other requirements to be useful in practice, such as having an acceptable variance.
#' In this work, we provide results on the operating characteristics of different estimators in optimal adaptive two-stage design
#' with normally distributed outcome. In optimal designs, design parameters such as the sample sizes, decision boundaries,
#' and adaptation rules are chosen as the result of an optimization process. The goal of the optimization process is to
#' maximize some metric of design quality, a typical example being the expected sample size required to fulfill certain
#' power requirements under a specified hypothesis. Although optimal adaptive designs have been a topic of recent research,
#' guidance on estimation in this novel kind of designs is still scarce. We compare classical and recently developed point
#' estimators as well as estimation methods for confidence intervals and p-values regarding various performance criteria such
#' as bias, variance, mean squared error, coverage, and consistency with test decisions.
#'
#' @docType package
#' @name adestr
#' @import adoptr methods
#' @importFrom stats dnorm pnorm qnorm dt pt qt dchisq pchisq qchisq integrate uniroot var
#' @importFrom cubature hcubature
## usethis namespace: start
#' @useDynLib adestr, .registration = TRUE
## usethis namespace: end
NULL

.adestr_options <- list(
  # Root finding inside estimators
  adestr_tol_roots = 1e-3,
  adestr_maxiter_roots = 1e3,
  # Integrals used inside estimators
  adestr_tol_inner = 5e-3,
  adestr_maxEval_inner = 1e4,
  adestr_absError_inner = 1e-5,
  # Integrals to evaluate estimators
  adestr_tol_outer = 5e-6,
  adestr_maxEval_outer = 3e5,
  adestr_absError_outer = 1e-8
)

#### Universal argument order for all functions: ####
# __Integration parameter
# x
# __User facing function arguments of generics:
# data, score, estimator, data_distribution, design, true_parameter, g0, g1,
# __Sufficient statistic parameters
# smean1, smean1T, smean1C, svar1, smean2, smean1T, smean1C, svar2, n1, n2, mu, sigma, two_armed,
# __Tuning parameters for various estimators
# d, p_boundary, alpha, beta, wc1f, wc1e, wc2, shiftc1f, shiftc1e, shiftc2
# __Other parameters
# x1, sample_space, smean1_region, combine_components
# __Parameters controlling quadrature and root finding last
# tol, maxEval, absError, exact, statistics
# __Parameters controlling plots
# facets
