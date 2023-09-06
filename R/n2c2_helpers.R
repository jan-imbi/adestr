n2_preimage <- function(design, sigma = 1, two_armed = FALSE, smean_scale = FALSE){
  design <- TwoStageDesignWithCache(design)
  zf <- design@c1f
  ze <- design@c1e
  n1 <- ceiling(n1(design, round = FALSE))
  se1 <- sigma * sqrt((1L + two_armed) / n1)
  x_candidates <- seq(design@c1f, design@c1e, length.out = 1e4)
  stepsize <- x_candidates[2L]-x_candidates[1L]
  n_candidates <- n2_extrapol(design, x_candidates)
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
            (ns[i] - 1L + sgn * .Machine$double.eps^.6) - n2_extrapol(design, x)
          },
          c(x_candidates[csum[i-1L]], x_candidates[csum[i-1L]] + stepsize),
          tol = .Machine$double.eps^.6
        )
      } else {
        root <- uniroot(
          function(x) {
            (ns[i] + sgn * .Machine$double.eps^.6) - n2_extrapol(design, x)
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
get_n2_coefficients <- function(design){
  h <- (design@c1e - design@c1f) / 2
  return(fastmonoH.FC_coefficients(
    h * design@x1_norm_pivots + (h + design@c1f),
    design@n2_pivots
  ))
}
get_c2_coefficients <- function(design){
  h <- (design@c1e - design@c1f) / 2
  return(fastmonoH.FC_coefficients(
    h * design@x1_norm_pivots + (h + design@c1f),
    design@c2_pivots
  ))
}


#' Calculate the second-stage sample size for a design with cached spline parameters
#'
#' Also extrapolates results for values outside of [c1f, c1e].
#'
#' @param design an object of class \code{\link{TwoStageDesignWithCache}}.
#' @param x1 first-stage test statistic
#'
#' @seealso \link[adoptr]{n2}
n2_extrapol <- function(design, x1) {
  if (length(design@n2_pivots > 1L)){
    fastmonoH.FC_evaluate(x1, design@n2_coefficients)
  } else {
    design@n2_pivots
  }
}
#' Calculate the second-stage critical value for a design with cached spline parameters
#'
#' Also extrapolates results for values outside of [c1f, c1e].
#'
#' @param design an object of class \code{\link{TwoStageDesignWithCache}}.
#' @param x1 first-stage test statistic
#'
#' @seealso \link[adoptr]{c2}
c2_extrapol <- function(design, x1) {
  fastmonoH.FC_evaluate(x1, design@c2_coefficients)
}

get_c2_extrapol_function <- function(design){
  h <- (design@c1e - design@c1f) / 2
  return(fastmonoH.FC_function(
    h * design@x1_norm_pivots + (h + design@c1f),
    design@c2_pivots
  ))
}
get_n2_extrapol_function <- function(design){
  if (length(design@n2_pivots)>1){
    h <- (design@c1e - design@c1f) / 2
    return(fastmonoH.FC_function(
      h * design@x1_norm_pivots + (h + design@c1f),
      design@n2_pivots
    ))
  } else{
    return(\(x) design@n2_pivots)
  }
}
