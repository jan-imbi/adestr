setClass(
  "TwoStageDesignWithCache",
  contains = "TwoStageDesign",
  slots = c(n2_coefficients = "list",
            c2_coefficients = "list")
)
TwoStageDesignWithCache <- function(design){
  new("TwoStageDesignWithCache",
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
}
forget_cache <- function(design){
  new("TwoStageDesign",
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
setMethod("print", signature("TwoStageDesignWithCache"), function(x, ...){
  print(forget_cache(x))
})


