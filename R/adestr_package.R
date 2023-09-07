#' adestr
#' @description Point estimates, confidence intervals, and p-values for optimal adaptive two-stage designs
#'
#' @details This package implements methods to evaluate the performance characteristics of
#' various point and interval estimators for optimal adaptive two-stage designs.
#' Specifically, this package is written to interface with trial designs created by the \code{\link{adoptr}} package
#' \insertCite{kunzmann2021adoptr,pilz2021optimal}{adestr}.
#' Apart from the a priori evaluation of performance characteristics, this package also allows for the
#' evaluation of the implemented estimators on real datasets, and it implements methods
#' to calculate p-values.
#'
#' @docType package
#' @name adestr
#' @import adoptr methods
#' @importFrom stats dnorm pnorm qnorm dt pt qt dchisq pchisq qchisq integrate uniroot var
#' @importFrom cubature hcubature
#' @importFrom Rdpack reprompt
#'
#' @references
#' \insertAllCited{}
## usethis namespace: start
#' @useDynLib adestr, .registration = TRUE
## usethis namespace: end
"_PACKAGE"

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
