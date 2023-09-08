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
setMethod("print", signature("TwoStageDesignWithCache"), function(x, ...){
  print(forget_cache(x))
})
