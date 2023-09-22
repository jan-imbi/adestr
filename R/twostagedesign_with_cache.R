### Remove some of this once adoptr is back on CRAN ###

#' Two-stage designs
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
TwoStageDesign <- function(n1, c1f, c1e, n2_pivots, c2_pivots, order = NULL) {
  # The original version of this is available from
  # https://github.com/kkmann/adoptr/blob/master/R/TwoStageDesign.R
  if (is.null(order)) {
    order <- length(c2_pivots)
  } else if (length(n2_pivots) != order) {
    n2_pivots <- rep(n2_pivots[1], order)
    c2_pivots <- rep(c2_pivots[1], order)
  }
  order <- as.integer(order)
  if (order < 2) stop("At least two nodes are necessary for integration!")
  j   <- 1:(order - 1)
  mu0 <- 2
  b   <- j / (4 * j^2 - 1)^0.5
  A   <- rep(0, order * order)
  A[(order + 1) * (j - 1) + 2] <- b
  A[(order + 1) * j] <- b
  dim(A) <- c(order, order)
  sd <- eigen(A, symmetric = TRUE)
  w <- rev(as.vector(sd$vectors[1, ]))
  w <- mu0 * w^2
  x <- rev(sd$values)
  rule <- data.frame(nodes = x, weights = w)
  tunable        <- logical(8) # initialize to all false
  tunable[1:5]   <- TRUE
  names(tunable) <- c("n1", "c1f", "c1e", "n2_pivots", "c2_pivots", "x1_norm_pivots", "weights", "tunable")
  new("TwoStageDesign", n1 = n1, c1f = c1f, c1e = c1e, n2_pivots = n2_pivots,
      c2_pivots = c2_pivots, x1_norm_pivots = rule$nodes, weights = rule$weights,
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
#'
#' @export
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
