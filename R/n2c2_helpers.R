n2_preimage <- function(design, sigma = 1, two_armed = FALSE, smean_scale = FALSE){
  zf <- design@c1f
  ze <- design@c1e
  n1 <- ceiling(n1(design, round = FALSE))
  se1 <- sigma * sqrt((1L + two_armed) / n1)
  x_candidates <- seq(design@c1f, design@c1e, length.out = 1e4)
  stepsize <- x_candidates[2L]-x_candidates[1L]
  n_candidates <- n2(design, x_candidates, round = FALSE)
  n_rle <- rle(ceiling(n_candidates))
  csum <- cumsum(n_rle$lengths)
  ns <- n_rle$values
  xs <- numeric(length(ns))
  for (i in seq_along(ns)){
    if (i == 1L){
      xs[i] <- zf
    } else{
      sgn <- sign(ns[i] - ns[i-1L])
      if (sgn>0){
        root <- uniroot(
          function(x) {
            (ns[i] - 1L + sgn * .Machine$double.eps^.6) - n2(design, x, round = FALSE)
          },
          c(x_candidates[csum[i-1L]], x_candidates[csum[i-1L]] + stepsize),
          tol = .Machine$double.eps^.6
        )
      } else {
        root <- uniroot(
          function(x) {
            (ns[i] + sgn * .Machine$double.eps^.6) - n2(design, x, round = FALSE)
          },
          c(x_candidates[csum[i-1L]], x_candidates[csum[i-1L]] + stepsize),
          tol = .Machine$double.eps^.6
        )
      }
      xs[i] <- root$root
    }
  }
  ret <- list()
  mult <- if (smean_scale) se1 else 1
  for (i in seq_along(xs)) {
    if (i < length(xs))
      ret[[i]] <- list(preimage = c(xs[i], xs[i+1]) * mult, n2 = ns[[i]])
    else
      ret[[i]] <- list(preimage = c(xs[i], ze) * mult, n2 = ns[[i]])
  }
  names(ret) <- as.character(ns)
  ret
}
cache_design_splines <- function(design, force = TRUE) {
  if (force | is.null(attr(design, "n2_cache")) | is.null("c2_cache")){
    attr(design, "n2_cache") <- get_fast_n2(design)
    attr(design, "c2_cache") <- get_fast_c2(design)
  }
  design
}
get_fast_c2 <- function(design){
  h <- (design@c1e - design@c1f) / 2
  return(fastmonoH.FC(
    h * design@x1_norm_pivots + (h + design@c1f),
    design@c2_pivots
  ))
}
get_fast_n2 <- function(design){
  if (length(design@n2_pivots)>1){
    h <- (design@c1e - design@c1f) / 2
    return(fastmonoH.FC(
      h * design@x1_norm_pivots + (h + design@c1f),
      design@n2_pivots
    ))
  } else{
    return(\(x) design@n2_pivots)
  }
}
n2_extrapol <- function(design, x1) {
  attr(design, "n2_cache")(x1)
}
c2_extrapol <- function(design, x1) {
  attr(design, "c2_cache")(x1)
}
# n2_extrapol <- function(design, x1) {
#   if (length(design@n2_pivots)>1){
#     h <- (design@c1e - design@c1f) / 2
#     return(stats::splinefun(
#       h * design@x1_norm_pivots + (h + design@c1f),
#       design@n2_pivots,
#       method = "monoH.FC"
#     )(x1))
#   } else{
#     return(design@n2_pivots)
#   }
# }
# c2_extrapol <- function(design, x1) {
#   h <- (design@c1e - design@c1f) / 2
#   return(stats::splinefun(
#     h * design@x1_norm_pivots + (h + design@c1f),
#     design@c2_pivots,
#     method = "monoH.FC"
#   )(x1))
# }
