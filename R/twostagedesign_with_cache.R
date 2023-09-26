### Remove some of this once adoptr is back on CRAN ###

#' Re-export of two-stage design class
#'
#' This is a re-export of the \code{TwoStageDesign} class from the
#' \code{adoptr} \insertCite{kunzmann2021adoptr}{adestr} package.
#'
#' This function is currently re-exported here to resolve CRAN conflicts.
#' For details, please refer to the paper \insertCite{kunzmann2021adoptr}{adestr}
#' and the original github repository \url{https://github.com/kkmann/adoptr}.
#'
#' @slot n1 (numeric) first-stage sample size.
#' @slot c1f (numeric) first-stage futility boundary.
#' @slot c1e (numeric) first-stage early efficacy boundary.
#' @slot n2_pivots (numeric) vector containing the values of the n2 spline function.
#' @slot c2_pivots (numeric) vector containing the values of the second-stage
#' rejection boundary spline c2
#' @slot x1_norm_pivots (numeric) vector containing the x-axis (z-sclae) points
#' for the n2 and c2 splines
#' @slot weights (numeric) vector containing integration weights
#' @slot tunable (logical) vector determining whether desing paramters are to be optimized
#'
#' @seealso The original implementation of the adoptr package by Kevin Kunzmann and
#' Maximilian Pilz is available at \url{https://github.com/kkmann/adoptr}.
#'
#' @exportClass TwoStageDesign
#' @importFrom pracma gaussLegendre
setClass("TwoStageDesign", representation(
  n1        = "numeric",
  c1f       = "numeric",
  c1e       = "numeric",
  n2_pivots = "numeric",
  c2_pivots = "numeric",
  x1_norm_pivots = "numeric",
  weights   = "numeric",
  tunable   = "logical"
))
TwoStageDesign <- function(n1, c1f, c1e, n2_pivots, c2_pivots, order = length(c2_pivots)) {
  if (order != length(c2_pivots)) {
    stop("order needs to be same length as c2_pivots")
  }
  if (length(n2_pivots) != length(c2_pivots) && length(n2_pivots)!=1) {
    stop("n2_pivots needs to be the same length as c2_pivots or of length 1.")
  }
  if (length(n2_pivots)==1L)
    n2_pivots <- rep(n2_pivots, order)
  glr <- gaussLegendre(order, -1, 1)
  tunable        <- c("n1" = TRUE, "c1f" = TRUE, "c1e" = TRUE, "n2_pivots" = TRUE, "c2_pivots" = TRUE,
                      "x1_norm_pivots" = FALSE, "weights" = FALSE, "tunable" = FALSE)
  new("TwoStageDesign", n1 = n1, c1f = c1f, c1e = c1e, n2_pivots = n2_pivots,
      c2_pivots = c2_pivots, x1_norm_pivots = glr$x, weights = glr$w,
      tunable = tunable)
}
setClass("DataDistribution", representation(
  two_armed = "logical")
)
setClass("Normal", contains = "DataDistribution")
#' Normally distributed data with known variance
#'
#' This function creates an object representing the distributional assumptions
#' of the data: normally distributed outcomes sample from a trial with
#' one or two arms (depending on the value of the parameter \code{two_armed}),
#' under the assumption of known variance.
#'
#' @param two_armed (logical) determines whether one or two-armed trials are assumed.
#'
#' @export
#' @returns an object of class \code{Normal}. This object encodes the distributional
#' assumptions of the data for usage in the functions
#' \code{\link{evaluate_estimator}} and \code{\link{analyze}}.
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
#'
Normal <- function(two_armed = TRUE) new("Normal", two_armed = two_armed)
setClass("Student", contains = "DataDistribution")

#' Normally distributed data with unknown variance
#'
#' This function creates an object representing the distributional assumptions
#' of the data: normally distributed outcomes sample from a trial with
#' one or two arms (depending on the value of the parameter \code{two_armed}),
#' under the assumption of known variance.
#'
#' @param two_armed (logical) determines whether one or two-armed trials are assumed.
#' @returns an object of class \code{Student}. This object encodes the distributional
#' assumptions of the data for usage in the functions
#' \code{\link{evaluate_estimator}} and \code{\link{analyze}}.
#'
#' @export
#' @examples
#' evaluate_estimator(
#'   score = MSE(),
#'   estimator = SampleMean(),
#'   data_distribution = Student(FALSE),
#'   design = get_example_design(),
#'   mu = c(0, 0.3, 0.6),
#'   sigma = 1,
#'   exact = FALSE
#' )
#'
Student <- function(two_armed = TRUE) new("Student", two_armed = two_armed)
n1 <- function(design, round = FALSE) if (round) round(design@n1) else design@n1
### end of remove ###

setClass(
  "TwoStageDesignWithCache",
  contains = "TwoStageDesign",
  slots = c(n2_coefficients = "list",
            c2_coefficients = "list")
)
#' TwoStageDesignWithCache constructor function
#'
#' Creates an object of class \code{TwoStageDesignWithCache}.
#' This object stores the precalculated spline paramters of the \code{n2}
#' and \code{c2} functions, which allows for quicker evaluation.
#'
#' @param design an object of class TwoStageDesign
#'
TwoStageDesignWithCache <- function(design) {
  if (attr(design, "class") == "TwoStageDesignWithCache") {
    return(design)
  } else {
    d <- new("TwoStageDesignWithCache",
      n1 = design@n1,
      c1f = design@c1f,
      c1e = design@c1e,
      n2_pivots = design@n2_pivots,
      c2_pivots = design@c2_pivots,
      x1_norm_pivots = design@x1_norm_pivots,
      weights = design@weights,
      tunable = design@tunable,
      n2_coefficients = get_n2_coefficients(design),
      c2_coefficients = get_c2_coefficients(design)
    )
    attr(d, "label") <- attr(design, "label")
    return(d)
  }
}
forget_cache <- function(design){
  label <- attr(design, "label")
  if (length(design@n2_pivots)==1) {
    d <- new("GroupSequentialDesign",
        n1 = design@n1,
        c1f = design@c1f,
        c1e = design@c1e,
        n2_pivots = design@n2_pivots,
        c2_pivots = design@c2_pivots,
        x1_norm_pivots = design@x1_norm_pivots,
        weights = design@weights,
        tunable = design@tunable
    )
  } else {
    d <- new("TwoStageDesign",
        n1 = design@n1,
        c1f = design@c1f,
        c1e = design@c1e,
        n2_pivots = design@n2_pivots,
        c2_pivots = design@c2_pivots,
        x1_norm_pivots = design@x1_norm_pivots,
        weights = design@weights,
        tunable = design@tunable
    )
  }
  attr(d, "label") <- label
  d
}
