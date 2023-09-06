# TODO: remove dots with this regex
# \.(?=[a-z|A-Z])

## Some helper functions.
.x1_x2_to_x <- function(x1, x2, n1, n2){
  (n1 * x1 + n2 * x2) / (n1 + n2)
}
.x_to_z<- function(x, n, mu0, sigma) {
  (x - mu0) * sqrt(n) / sigma
}
.x1_x2_to_z <- function(x1, x2, n1, n2, mu0, sigma) {
  (.x1_x2_to_x(x1 = x1, x2 = x2, n1 = n1, n2 = n2) - mu0) * sqrt(n1 + n2) / sigma
}
.z_to_x <- function(z, n, mu0, sigma){
  z * sigma / sqrt(n) + mu0
}
.z1_z2_to_x <- function(z1, z2, n1, n2, mu0, sigma){
  x1 <- .z_to_x(z = z1, n = n1, mu0 = mu0, sigma = sigma)
  x2 <- .z_to_x(z = z2, n = n2, mu0 = mu0, sigma = sigma)
  (n1 * x1 + n2 * x2) / (n1 + n2)
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

## The densities for integration.
.f1 <- function(z1, n1, mu, mu0, sigma) {
  dnorm(z1 + (mu0 - mu) * sqrt(n1)/sigma)
}
.f2 <- function(z1, z2, n1, n2, mu, mu0, sigma){
  .f1(z1 = z1, n1 = n1, mu = mu, mu0 = mu0, sigma = sigma) * dnorm(z2 + (mu0 - mu) * sqrt(n2)/sigma)
}
.f <- function(z1, z2, mu, mu0, sigma, design) {
  n1 <- design@n1
  c1f <- design@c1f
  c1e <- design@c1e
  if ((z1 < c1f || z1 > c1e) && (z2 != -Inf)) {
    stop("z1 suggests early stopping but z2 is not -Inf.")
  }
  if ((z1 >= c1f  && z1 <= c1e) && (z2 == -Inf)) {
    stop("z1 suggests continuation but z2 is -Inf.")
  }
  ret <-
    if (z1 < c1f || z1 > c1e) {
      .f1(z1 = z1, n1 = n1, mu = mu, mu0 = mu0, sigma = sigma)
    } else {
      n2 <- n2_extrapol(design, z1)
      .f2(z1 = z1, z2 = z2, n1 = n1, n2 = n2, mu = mu, mu0 = mu0, sigma = sigma)
    }
  ret
}

## Point estimators
### Maximum likelihood estimator (MLE) and functions to evaluate its distribution,
### expectation and variance
.mle <- function(x1, x2, n1, n2) {
  (n1 * x1 + n2 * x2) / (n1 + n2)
}
.mle_cdf <- function(y, mu, mu0, sigma, design) {
  n1 <- design@n1
  yz <- .x_to_z(x = y, n = n1, mu0 = mu0, sigma = sigma)
  c1f <- design@c1f
  c1e <- design@c1e
  early_stopping_int <-
    if (yz < design@c1f) {
      hcubature(
        .f1,
        vectorInterface = TRUE,
        lowerLimit = -Inf,
        upperLimit = yz,
        n1 = n1, mu = mu,  mu0 = mu0, sigma = sigma
      )$integral
    } else if (yz > design@c1e) {
      lower_int <- hcubature(
        .f1,
        vectorInterface = TRUE,
        lowerLimit = -Inf,
        upperLimit = c1f,
        n1 = n1, mu = mu,  mu0 = mu0, sigma = sigma
      )$integral
        upper_int <- hcubature(
        .f1,
        vectorInterface = TRUE,
        lowerLimit = c1e,
        upperLimit = yz,
        n1 = n1, mu = mu,  mu0 = mu0, sigma = sigma
      )$integral
      lower_int + upper_int
    } else {
      hcubature(
        .f1,
        vectorInterface = TRUE,
        lowerLimit = -Inf,
        upperLimit = c1f,
        n1 = n1, mu = mu,  mu0 = mu0, sigma = sigma
      )$integral
    }
  # Cannot use 2d integral because boundaries of integral are function of z1
  continuation_int <- hcubature(
    \(z1) {
       n2 <- n2_extrapol(design, z1)
       x1 <- .z_to_x(z = z1, n = n1, mu0 = mu0, sigma = sigma)
       .f1(z1 = z1, n1 = n1, mu = mu, mu0 = mu0, sigma = sigma) *
       hcubature(
         \(z2, n2) .f1(z1 = z2[1,,drop=FALSE], n1 = n2, mu = mu, mu0 = mu0, sigma = sigma),
         lowerLimit = -Inf,
         upperLimit = (((n1 + n2) * y - n1 * x1)/n2 - mu0) * sqrt(n2)/sigma,
         vectorInterface = TRUE,
         n2 = n2
       )$integral
    },
    vectorInterface = FALSE,
    lowerLimit = c(design@c1f),
    upperLimit = c(design@c1e)
  )$integral
  early_stopping_int + continuation_int
}
.mle_pdf <- function(y, mu, mu0, sigma, design) {
  n1 <- design@n1
  int <- hcubature(
    \(x, n1, mu, mu0, sigma, design) {
      n2 <- n2_extrapol(design, x[1,,drop=FALSE])
      z2 <- ( ((n1 + n2) * y - n1 * .z_to_x(z = x[1,,drop=FALSE], n = n1, mu0 = mu0, sigma = sigma)) /n2 - mu0) * sqrt(n2)/sigma
      (n1 + n2)/n2 * sqrt(n2)/sigma * .f2(z1 = x[1,,drop=FALSE], z2 = z2, n1 = n1, n2 = n2, mu = mu, mu0 = mu0, sigma = sigma)
    },
    vectorInterface = TRUE,
    lowerLimit = design@c1f,
    upperLimit = design@c1e,
    n1 = n1, mu = mu, mu0 = mu0, sigma = sigma, design = design
  )$integral
  yz1 <- (y - mu0) * sqrt(n1)/sigma
  second_part <-
  if (yz1 < design@c1f || yz1 > design@c1e)
    sqrt(n1)/sigma * .f1(z1 = yz1, n1 = n1, mu = mu, mu0 = mu0, sigma = sigma)
  else
    0
  int + second_part
}
.mle_expectation <- function(mu, mu0, sigma, design){
  n1 <- design@n1
  futility <- hcubature(
    f = \(x, n1, mu, mu0, sigma){
      .z_to_x(z = x[1,,drop=FALSE], n = n1, mu0 = mu0, sigma = sigma) *
        .f1(z1 = x[1,,drop=FALSE], n1 = n1, mu = mu, mu0 = mu0, sigma = sigma)
    },
    lowerLimit = c(-Inf),
    upperLimit = c(design@c1f),
    vectorInterface = TRUE,
    n1 = n1, mu = mu, mu0 = mu0, sigma = sigma
  )$integral
  efficacy <- hcubature(
    f = \(x, n1, mu, mu0, sigma){
      .z_to_x(z = x[1,,drop=FALSE], n = n1, mu0 = mu0, sigma = sigma) *
        .f1(z1 = x[1,,drop=FALSE], n1 = n1, mu = mu, mu0 = mu0, sigma = sigma)
    },
    lowerLimit = c(design@c1e),
    upperLimit = c(Inf),
    vectorInterface = TRUE,
    n1 = n1, mu = mu, mu0 = mu0, sigma = sigma
  )$integral
  continuation <- hcubature(
    f = \(x, n1, mu, mu0, sigma, design){
      n2 <- n2_extrapol(design, x[1,,drop=FALSE])
      .z1_z2_to_x(z1 = x[1, , drop = FALSE], z2 = x[2, , drop = FALSE], n1 = n1, n2 = n2, mu0 = mu0, sigma = sigma ) *
        .f2(z1 = x[1, , drop = FALSE], z2 = x[2, , drop = FALSE], n1 = n1, n2 = n2, mu = mu, mu0 = mu0, sigma = sigma)
    },
    lowerLimit = c(design@c1f,-Inf),
    upperLimit = c(design@c1e, Inf),
    vectorInterface = TRUE,
    n1 = n1, mu = mu, mu0 = mu0, sigma = sigma, design = design
  )$integral
  efficacy + futility + continuation
}

### Fixed and adaptivly weighted sample means.
.fixed_weighted_mle <- function(x1, x2, w,...){
  w * x1 + (1-w) * x2
}
.adaptively_weighted_mle <- function(x1, x2, n1, n2, sigma, w){
  tau <- (w * sqrt(n1) / sigma) / (w * sqrt(n1) / sigma + sqrt(1 - w^2) * sqrt(n2) / sigma)
  tau * x1 + (1 - tau) * x2
}
.adaptively_weighted_mle_expectation <- function(w, mu, mu0, sigma, design) {
  n1 <- design@n1
  futility <- hcubature(
    f = \(x, n1, mu, mu0, sigma){
      .z_to_x(z = x[1,,drop=FALSE], n = n1, mu0 = mu0, sigma = sigma) *
        .f1(z1 = x[1,,drop=FALSE], n1 = n1, mu = mu, mu0 = mu0, sigma = sigma)
    },
    lowerLimit = c(-Inf),
    upperLimit = c(design@c1f),
    vectorInterface = TRUE,
    n1 = n1, mu = mu, mu0 = mu0, sigma = sigma
  )$integral
  efficacy <- hcubature(
    f = \(x, n1, mu, mu0, sigma){
      .z_to_x(z = x[1,,drop=FALSE], n = n1, mu0 = mu0, sigma = sigma) *
        .f1(z1 = x[1,,drop=FALSE], n1 = n1, mu = mu, mu0 = mu0, sigma = sigma)
    },
    lowerLimit = c(design@c1e),
    upperLimit = c(Inf),
    vectorInterface = TRUE,
    n1 = n1, mu = mu, mu0 = mu0, sigma = sigma
  )$integral
  continuation <- hcubature(
    f = \(x, n1, mu, mu0, sigma, design){
      n2 <- n2_extrapol(design, x[1,,drop=FALSE])
      x1 <- .z_to_x(z = x[1,,drop=FALSE], n = n1, mu0 = mu0, sigma = sigma)
      x2 <- .z_to_x(z = x[2,,drop=FALSE], n = n2, mu0 = mu0, sigma = sigma)
      .adaptively_weighted_mle(x1 = x1, x2 = x2, n1 = n1, n2 = n2, sigma = sigma, w = w) *
        .f2(z1 = x[1, , drop = FALSE], z2 = x[2, , drop = FALSE], n1 = n1, n2 = n2, mu = mu, mu0 = mu0, sigma = sigma)
    },
    lowerLimit = c(design@c1f,-Inf),
    upperLimit = c(design@c1e, Inf),
    vectorInterface = TRUE,
    n1 = n1, mu = mu, mu0 = mu0, sigma = sigma, design = design
  )$integral
  efficacy + futility + continuation
}
.adaptively_weighted_mle_variance <- function(w, mu, mu0, sigma, design) {
  E <- .adaptively_weighted_mle_expectation(w, mu, mu0, sigma, design)
  n1 <- design@n1
  futility <- hcubature(
    f = \(x, n1, mu, mu0, sigma){
      (.z_to_x(z = x[1,,drop=FALSE], n = n1, mu0 = mu0, sigma = sigma) - E)^2 *
        .f1(z1 = x[1,,drop=FALSE], n1 = n1, mu = mu, mu0 = mu0, sigma = sigma)
    },
    lowerLimit = c(-Inf),
    upperLimit = c(design@c1f),
    vectorInterface = TRUE,
    n1 = n1, mu = mu, mu0 = mu0, sigma = sigma
  )$integral
  efficacy <- hcubature(
    f = \(x, n1, mu, mu0, sigma){
      (.z_to_x(z = x[1,,drop=FALSE], n = n1, mu0 = mu0, sigma = sigma) - E)^2 *
        .f1(z1 = x[1,,drop=FALSE], n1 = n1, mu = mu, mu0 = mu0, sigma = sigma)
    },
    lowerLimit = c(design@c1e),
    upperLimit = c(Inf),
    vectorInterface = TRUE,
    n1 = n1, mu = mu, mu0 = mu0, sigma = sigma
  )$integral
  continuation <- hcubature(
    f = \(x, n1, mu, mu0, sigma, design){
      n2 <- n2_extrapol(design, x[1,,drop=FALSE])
      x1 <- .z_to_x(x[1,,drop=FALSE], n = n1, mu0 = mu0, sigma = sigma)
      x2 <- .z_to_x(x[2,,drop=FALSE], n = n2, mu0 = mu0, sigma = sigma)
      (.adaptively_weighted_mle(x1 = x1, x2 = x2, n1 = n1, n2 = n2, sigma = sigma, w = w) - E)^2 *
        .f2(z1 = x[1,,drop = FALSE], z2 = x[2,,drop = FALSE], n1 = n1, n2 = n2, mu = mu, mu0 = mu0, sigma = sigma)
    },
    lowerLimit = c(design@c1f,-Inf),
    upperLimit = c(design@c1e, Inf),
    vectorInterface = TRUE,
    n1 = n1, mu = mu, mu0 = mu0, sigma = sigma, design = design
  )$integral
  efficacy + futility + continuation
}
.determine_min_variance_weight <- function(mu, mu0, sigma, design, weight_range = c(1e-10, 1-1e-10)) {
  n1 <- design@n1
  minw <- optimize(
    .adaptively_weighted_mle_variance,
    interval = weight_range,
    maximum = FALSE,
    mu = mu,
    mu0 = mu0,
    sigma = sigma,
    design = design
  )$minimum
  max_n2 <- n2_extrapol(design, design@c1f)
  min_n2 <- n2_extrapol(design, design@c1e)
  max_tau <- (minw * sqrt(design@n1) / sigma) / (minw * sqrt(design@n1) / sigma + sqrt(1 - minw^2) * sqrt(min_n2) / sigma)
  min_tau <- (minw * sqrt(design@n1) / sigma) / (minw * sqrt(design@n1) / sigma + sqrt(1 - minw^2) * sqrt(max_n2) / sigma)
  message(sprintf("Minimum weight w = %.3f, which means tau ranges from %.3f to %.3f.", minw, min_tau, max_tau))
  minw
}

### (Pseudo) Rao-Blackwell estimators.
.pseudo_rao_blackwell <- function(x1, x2, mu0, sigma, design) {
  n1 <- design@n1
  c1f <- design@c1f
  c1e <- design@c1e
  z1 <- .x_to_z(x = x1, n = n1, mu0 = mu0, sigma = sigma)
  if ((z1 < c1f || z1 > c1e) && (x2 != -Inf)) {
    stop("z1 suggests early stopping but x2 is not -Inf.")
  }
  if ((z1 >= c1f && z1 <= c1e) && (x2 == -Inf)) {
    stop("z1 suggests continuation but x2 is -Inf.")
  }
  ret <-
    if (z1 < c1f || z1 > c1e) {
      x1
    } else {
      n2 <- n2_extrapol(design, z1)
      x <- .x1_x2_to_x(x1 = x1, x2 = x2, n1 = n1, n2 = n2)
      numerator <-
        hcubature(
          f = \(z1_) {
            x1_ <- .z_to_x(z = z1_[1,,drop=FALSE], n = n1, mu0 = mu0, sigma = sigma)
            x1_ * .f2(z1 = z1_[1,,drop=FALSE], z2 = (((n1 + n2) * x - n1 * x1_) / n2 - mu0) * sqrt(n2) / sigma,
                      n1 = n1,  n2 = n2, mu = mu0, mu0 = mu0, sigma = sigma)
          },
          lowerLimit = c1f,
          upperLimit = c1e,
          vectorInterface = TRUE
        )$integral
      denominator <-
        hcubature(
          f = \(z1_) {
            x1_ <- .z_to_x(z = z1_[1,,drop=FALSE], n = n1, mu0 = mu0, sigma = sigma)
            .f2(z1 = z1_[1,,drop=FALSE], z2 = (((n1 + n2) * x - n1 * x1_) / n2 - mu0) * sqrt(n2) / sigma,
                n1 = n1,  n2 = n2, mu = mu0, mu0 = mu0, sigma = sigma)
          },
          lowerLimit = c1f,
          upperLimit = c1e,
          vectorInterface = TRUE
        )$integral
      numerator / denominator
    }
  ret
}
.rao_blackwell <- function(x1, x2, mu0, sigma, design) {
  n1 <- design@n1
  c1f <- design@c1f
  c1e <- design@c1e
  z1 <- .x_to_z(x = x1, n = n1, mu0 = mu0, sigma = sigma)
  if ((z1 < c1f || z1 > c1e) && (x2 != -Inf)) {
    stop("z1 suggests early stopping but x2 is not -Inf.")
  }
  if ((z1 >= c1f && z1 <= c1e) && (x2 == -Inf)) {
    stop("z1 suggests continuation but x2 is -Inf.")
  }
  ret <-
    if (z1 < c1f || z1 > c1e) {
      x1
    } else {
      # The difference between the pseudo Rao-Blackwell and the Rao-Blackwell
      # estimator is that the (normal) Rao-Blackwell estimator uses the actual
      # n2 preimage.
      n2 <- ceiling(n2_extrapol(design, z1))
      lower_z1 <- uniroot(\(x) n2_extrapol(design, x) - n2 - -2, c(c1f, c1e))$root
      upper_z1 <- uniroot(\(x) n2_extrapol(design, x) - n2 - -1, c(c1f, c1e))$root
      x <- .x1_x2_to_x(x1 = x1, x2 = x2, n1 = n1, n2 = n2)
      numerator <-
        hcubature(
          f = \(z1_) {
            x1_ <- .z_to_x(z = z1_[1,,drop=FALSE], n = n1, mu0 = mu0, sigma = sigma)
            x1_ * .f2(z1 = z1_[1,,drop=FALSE], z2 = (((n1 + n2) * x - n1 * x1_) / n2 - mu0) * sqrt(n2) / sigma,
                      n1 = n1,  n2 = n2, mu = mu0, mu0 = mu0, sigma = sigma)
          },
          lowerLimit = c1f,
          upperLimit = c1e,
          vectorInterface = TRUE
        )$integral
      denominator <-
        hcubature(
          f = \(z1_) {
            x1_ <- .z_to_x(z = z1_[1,,drop=FALSE], n = n1, mu0 = mu0, sigma = sigma)
            .f2(z1 = z1_[1,,drop=FALSE], z2 = (((n1 + n2) * x - n1 * x1_) / n2 - mu0) * sqrt(n2) / sigma,
                n1 = n1,  n2 = n2, mu = mu0, mu0 = mu0, sigma = sigma)
          },
          lowerLimit = lower_z1,
          upperLimit = upper_z1,
          vectorInterface = TRUE
        )$integral
      numerator / denominator
    }
  ret
}

### Bias reduced estimator.
.mle_expectation_derivative <- function(mu, mu0, sigma, design){
  n1 <- design@n1
  # (Compare .mle_expectation_derivative to .mle_expectation)
  futility <- hcubature(
    f = \(x, n1, mu, mu0, sigma){
      .z_to_x(z = x[1,,drop=FALSE], n = n1, mu0 = mu0, sigma = sigma) *
        (x[1,,drop=FALSE] + (mu0 - mu) * sqrt(n1)/sigma) * sqrt(n1)/sigma *
        .f1(z1 = x[1,,drop=FALSE], n1 = n1, mu = mu, mu0 = mu0, sigma = sigma)
    },
    lowerLimit = c(-Inf),
    upperLimit = c(design@c1f),
    vectorInterface = TRUE,
    n1 = n1, mu = mu, mu0 = mu0, sigma = sigma
  )$integral
  efficacy <- hcubature(
    f = \(x, n1, mu, mu0, sigma){
      .z_to_x(z = x[1,,drop=FALSE], n = n1, mu0 = mu0, sigma = sigma) *
        (x[1,,drop=FALSE] + (mu0 - mu) * sqrt(n1)/sigma) * sqrt(n1)/sigma *
        .f1(z1 = x[1,,drop=FALSE], n1 = n1, mu = mu, mu0 = mu0, sigma = sigma)
    },
    lowerLimit = c(design@c1e),
    upperLimit = c(Inf),
    vectorInterface = TRUE,
    n1 = n1, mu = mu, mu0 = mu0, sigma = sigma
  )$integral
  continuation <- hcubature(
    f = \(x, n1, mu, mu0, sigma, design){
      n2 <- n2_extrapol(design, x[1,,drop=FALSE])
      .z1_z2_to_x(z1 = x[1, , drop = FALSE], z2 = x[2, , drop = FALSE], n1 = n1, n2 = n2, mu0 = mu0, sigma = sigma) *
        ((x[1,,drop=FALSE] + (mu0 - mu) * sqrt(n1)/sigma) * sqrt(n1)/sigma +
           (x[2,,drop=FALSE] + (mu0 - mu) * sqrt(n2)/sigma) * sqrt(n2)/sigma) *
        .f2(z1 = x[1, , drop = FALSE], z2 = x[2, , drop = FALSE], n1 = n1, n2 = n2, mu = mu, mu0 = mu0, sigma = sigma)
    },
    lowerLimit = c(design@c1f,-Inf),
    upperLimit = c(design@c1e, Inf),
    vectorInterface = TRUE,
    n1 = n1, mu = mu, mu0 = mu0, sigma = sigma, design = design
  )$integral
  efficacy + futility + continuation
}
.bias_reduced <- function(mle, mu, mu0, sigma, iterations, design){
  estimate <- mle
  for (i in seq_len(iterations)) {
    expectation_mle <- .mle_expectation(mu = estimate, mu0 = mu0, sigma = sigma, design = design)
    derivative_mle <- .mle_expectation_derivative(mu = estimate, mu0 = mu0, sigma = sigma, design = design)
    estimate <- estimate +  ((mle - estimate) - (expectation_mle - mle)) / (1 + derivative_mle)
  }
  estimate
}

## Ordering based methods.
###  Sample space order inducing functions and helper functions for y1, A and B.
.rho_ml1 <- function(x1){
  x1
}
.rho_ml2 <- function(x1, x2, n1, n2){
  .x1_x2_to_x(x1 = x1, x2 = x2, n1 = n1, n2 = n2)
}
.rho_ml_y <- function(rho_ml){
  rho_ml
}
.rho_ml_B <- function(rho_ml, x1, n1, n2, mu){
  (rho_ml * (n1 + n2) - n1 * x1) / n2
}
.rho_lr1 <- function(x1, n1, mu) {
  (x1 - mu) * sqrt(n1)
}
.rho_lr2 <- function(x1, x2, n1, n2, mu) {
  (.x1_x2_to_x(x1 = x1, x2 = x2, n1 = n1, n2 = n2) - mu) * sqrt(n1 + n2)
}
.rho_lr_y <- function(rho_lr, n1, mu){
  (rho_lr + mu * sqrt(n1))/sqrt(n1)
}
.rho_lr_B <- function(rho_lr, x1, n1, n2, mu){
  (rho_lr * sqrt(n1 + n2) - n1 * x1 + mu * (n1 + n2))/n2
}
.rho_st1 <- function(x1, n1, mu) {
  n1 * x1 - n1 * mu
}
.rho_st2 <- function(x1, x2, n1, n2, mu) {
  n1 * x1 + n2 * x2 - (n1 + n2) * mu
}
.rho_st_y <- function(rho_st, n1, mu){
  (rho_st + mu*n1)/n1
}
.rho_st_B <- function(rho_st, x1, n1, n2, mu){
  (rho_st - n1 * x1 + mu * (n1 + n2))/n2
}
.rho_np1 <- function(x1, n1, mu0, mu1) {
  (2*x1*(mu0 - mu1) - mu0 - mu1)*n1
}
.rho_np2 <- function(x1, x2, n1, n2, mu0, mu1) {
  (2*.x1_x2_to_x(x1 = x1, x2 = x2, n1 = n1, n2 = n2)*(mu0 - mu1) - mu0 - mu1)*(n1 + n2)
}
.rho_np_y <- function(rho_np, n1, mu0, mu1){
  (rho_np + (mu0 + mu1)*n1) / (2*n1)
}
.rho_np_B <- function(rho_np, x1, n1, n2, mu0, mu1){
  (rho_np - 2 * n1 * x1 + (mu0 + mu1)*(n1 + n2)) / (2*n1)
}
.rho_swcf1 <- function(x1, n1, mu0, sigma){
  pnorm(.x_to_z(x = x1, n = n1, mu0 = mu0, sigma = sigma))
}
.rho_swcf2 <- function(x1, x2, n1, n2, mu, sigma, design){
  1 + pnorm(.swcf(x1 = x1, x2 = x2, n1 = n1, n2 = n2, mu = mu, sigma = sigma, design = design))
}
.rho_swcf3 <- function(x1, n1, mu0, sigma){
  2 + pnorm(.x_to_z(x = x1, n = n1, mu0 = mu0, sigma = sigma))
}
.swcf <- function(x1, x2, n1, n2, mu, sigma, design){
  .x_to_z(x = x2, n = n2, mu0 = mu, sigma = sigma) - c2_extrapol(design, .x_to_z(x = x1, n = n1, mu0 = mu, sigma = sigma))
}
.rho_swcf_B <- function(cf_value, x1, n1, n2, mu, sigma, design){
  sigma / sqrt(n2) * (cf_value + c2_extrapol(design, (x1 - mu)*sqrt(n1)/sigma)) + mu
}
.calculate_A <- function(y, n1, mu, sigma, design){
  c1f <- design@c1f
  c1e <- design@c1e
  ret <-
    if (y < c1f * sigma / sqrt(n1)) {
      pnorm(c1f - mu * sqrt(n1)/sigma) -
        pnorm(y * sqrt(n1) / sigma - mu * sqrt(n1)/sigma) +
        1 - pnorm(c1e - mu*sqrt(n1)/sigma)
    } else if (y > c1e * sigma / sqrt(n1)) {
      1 - pnorm(y * sqrt(n1) / sigma - mu * sqrt(n1)/sigma)
    } else {
      1 -  pnorm(c1e - mu*sqrt(n1)/sigma)
    }
  ret
}

### MLE ordering
#### P-value (MLE ordering).
.p_ml <- function(x1, x2, mu, mu0, sigma, design){
  n1 <- design@n1
  c1f <- design@c1f
  c1e <- design@c1e
  z1 <- .x_to_z(x = x1, n = n1, mu0 = mu0, sigma = sigma)
  if ((z1 < c1f || z1 > c1e) && (x2 != -Inf)) {
    stop("z1 suggests early stopping but x2 is not -Inf.")
  }
  if ((z1 >= c1f  && z1 <= c1e) && (x2 == -Inf)) {
    stop("z1 suggests continuation but x2 is -Inf.")
  }
  rho <-
    if (z1 < c1f || z1 > c1e) {
      .rho_ml1(x1 = x1)
    } else {
      n2 <- n2_extrapol(design, z1)
      .rho_ml2(x1 = x1, x2 = x2, n1 = n1, n2 = n2)
    }
  y <- .rho_ml_y(rho_ml = rho)
  A <- .calculate_A(y = y, n1 = n1, mu = mu, sigma = sigma, design = design)
  int <- hcubature(
    f = \(z_) {
      n2_ <- n2_extrapol(design, z_)
      x1_ <- .z_to_x(z_, n = n1, mu0 = mu0, sigma = sigma)
      (1 - pnorm((.rho_ml_B(rho_ml = rho, x1 = x1_, n1 = n1, n2 = n2_, mu = mu) - mu) * sqrt(n2_) / sigma)) *
        .f1(z1 = z_, n1 = n1, mu = mu, mu0 = mu0, sigma = sigma)
    },
    lowerLimit = c1f,
    upperLimit = c1e,
    vectorInterface = TRUE
  )$integral
  int + A
}
.rho_ml_root <- function(gamma, x1, x2, mu0, sigma, design){
  uniroot(
    f = \(x) .p_ml(x1 = x1, x2 = x2, mu = x, mu0 = mu0, sigma = sigma, design = design) - gamma,
    interval = qnorm(c(.025, .975), mean = mu0, sd = sigma / sqrt(design@n1)),
    extendInt = "yes"
  )$root
}
#### Median unbiased point estimator (MLE ordering).
.median_unbiased_ml <- function(x1, x2, mu0, sigma, design){
  .rho_ml_root(gamma = 0.5, x1 = x1, x2 = x2, mu0 = mu0, sigma = sigma, design = design)
}
#### Confidence interval (MLE ordering).
.confidence_interval_ml <- function(alpha, x1, x2, mu0, sigma, design){
  c(.rho_ml_root(gamma = alpha/2, x1 = x1, x2 = x2, mu0 = mu0, sigma = sigma, design = design),
    .rho_ml_root(gamma = 1-alpha/2, x1 = x1, x2 = x2, mu0 = mu0, sigma = sigma, design = design))
}

### Likelihood ratio ordering.
#### P-value (LR ordering).
.p_lr <- function(x1, x2, mu, mu0, sigma, design){
  n1 <- design@n1
  c1f <- design@c1f
  c1e <- design@c1e
  z1 <- .x_to_z(x = x1, n = n1, mu0 = mu0, sigma = sigma)
  if ((z1 < c1f || z1 > c1e) && (x2 != -Inf)) {
    stop("z1 suggests early stopping but x2 is not -Inf.")
  }
  if ((z1 >= c1f  && z1 <= c1e) && (x2 == -Inf)) {
    stop("z1 suggests continuation but x2 is -Inf.")
  }
  rho <-
    if (z1 < c1f || z1 > c1e) {
      .rho_lr1(x1 = x1, n1 = n1, mu = mu)
    } else {
      n2 <- n2_extrapol(design, z1)
      .rho_lr2(x1 = x1, x2 = x2, n1 = n1, n2 = n2, mu = mu)
    }
  y <- .rho_lr_y(rho_lr = rho, n1 = n1, mu = mu)
  A <- .calculate_A(y = y, n1 = n1, mu = mu, sigma = sigma, design = design)
  int <- hcubature(
    f = \(z_) {
      n2_ <- n2_extrapol(design, z_)
      x1_ <- .z_to_x(z_, n = n1, mu0 = mu0, sigma = sigma)
      (1 - pnorm((.rho_lr_B(rho_lr = rho, x1 = x1_, n1 = n1, n2 = n2_, mu = mu) - mu) * sqrt(n2_) / sigma)) *
        .f1(z1 = z_, n1 = n1, mu = mu, mu0 = mu0, sigma = sigma)
    },
    lowerLimit = c1f,
    upperLimit = c1e,
    vectorInterface = TRUE
  )$integral
  int + A
}
.rho_lr_root <- function(gamma, x1, x2, mu0, sigma, design){
  uniroot(
    f = \(x) .p_lr(x1 = x1, x2 = x2, mu = x, mu0 = mu0, sigma = sigma, design = design) - gamma,
    interval = qnorm(c(.025, .975), mean = mu0, sd = sigma / sqrt(design@n1)),
    extendInt = "yes"
  )$root
}
#### Median unbiased estimator (LR ordering).
.median_unbiased_lr <- function(x1, x2, mu0, sigma, design){
  .rho_lr_root(gamma = 0.5, x1 = x1, x2 = x2, mu0 = mu0, sigma = sigma, design = design)
}
#### Confidence interval (LR ordering).
.confidence_interval_lr <- function(alpha, x1, x2, mu0, sigma, design){
  c(.rho_lr_root(gamma = alpha/2, x1 = x1, x2 = x2, mu0 = mu0, sigma = sigma, design = design),
    .rho_lr_root(gamma = 1-alpha/2, x1 = x1, x2 = x2, mu0 = mu0, sigma = sigma, design = design))
}

### Score test ordering.
#### P-value (ST ordering).
.p_st <- function(x1, x2, mu, mu0, sigma, design){
  n1 <- design@n1
  c1f <- design@c1f
  c1e <- design@c1e
  z1 <- .x_to_z(x = x1, n = n1, mu0 = mu0, sigma = sigma)
  if ((z1 < c1f || z1 > c1e) && (x2 != -Inf)) {
    stop("z1 suggests early stopping but x2 is not -Inf.")
  }
  if ((z1 >= c1f  && z1 <= c1e) && (x2 == -Inf)) {
    stop("z1 suggests continuation but x2 is -Inf.")
  }
  rho <-
    if (z1 < c1f || z1 > c1e) {
      .rho_st1(x1 = x1, n1 = n1, mu = mu)
    } else {
      n2 <- n2_extrapol(design, z1)
      .rho_st2(x1 = x1, x2 = x2, n1 = n1, n2 = n2, mu = mu)
    }
  y <- .rho_st_y(rho_st = rho)
  A <- .calculate_A(y = y, n1 = n1, mu = mu, sigma = sigma, design = design)
  int <- hcubature(
    f = \(z_) {
      n2_ <- n2_extrapol(design, z_)
      x1_ <- .z_to_x(z_, n = n1, mu0 = mu0, sigma = sigma)
      (1 - pnorm((.rho_st_B(rho_st = rho, x1 = x1_, n1 = n1, n2 = n2_, mu = mu) - mu) * sqrt(n2_) / sigma)) *
        .f1(z1 = z_, n1 = n1, mu = mu, mu0 = mu0, sigma = sigma)
    },
    lowerLimit = c1f,
    upperLimit = c1e,
    vectorInterface = TRUE
  )$integral
  int + A
}
.rho_st_root <- function(gamma, x1, x2, mu0, sigma, design){
  uniroot(
    f = \(x) .p_st(x1 = x1, x2 = x2, mu = x, mu0 = mu0, sigma = sigma, design = design) - gamma,
    interval = qnorm(c(.025, .975), mean = mu0, sd = sigma / sqrt(design@n1)),
    extendInt = "yes"
  )$root
}
#### Median unbiased point estimator (ST ordering)
.median_unbiased_st <- function(x1, x2, mu0, sigma, design){
  .rho_st_root(gamma = 0.5, x1 = x1, x2 = x2, mu0 = mu0, sigma = sigma, design = design)
}
#### Confidence interval (ST ordering)
.confidence_interval_st <- function(alpha, x1, x2, mu0, sigma, design){
  c(.rho_st_root(gamma = alpha/2, x1 = x1, x2 = x2, mu0 = mu0, sigma = sigma, design = design),
    .rho_st_root(gamma = 1-alpha/2, x1 = x1, x2 = x2, mu0 = mu0, sigma = sigma, design = design))
}

### Neyman-Pearson ordering
#### P-value (NP ordering).
.p_np <- function(x1, x2, mu, mu0, mu1, sigma, design){
  # (In reality, this only makes sense for mu = mu0. But for the sake of completeness,
  # those parameters are separated.)
  n1 <- design@n1
  c1f <- design@c1f
  c1e <- design@c1e
  z1 <- .x_to_z(x = x1, n = n1, mu0 = mu0, sigma = sigma)
  if ((z1 < c1f || z1 > c1e) && (x2 != -Inf)) {
    stop("z1 suggests early stopping but x2 is not -Inf.")
  }
  if ((z1 >= c1f  && z1 <= c1e) && (x2 == -Inf)) {
    stop("z1 suggests continuation but x2 is -Inf.")
  }
  rho <-
    if (z1 < c1f || z1 > c1e) {
      .rho_np1(x1 = x1, n1 = n1, mu0 = mu0, mu1 = mu1)
    } else {
      n2 <- n2_extrapol(design, z1)
      .rho_np2(x1 = x1, x2 = x2, n1 = n1, n2 = n2, mu0 = mu0, mu1 = mu1)
    }
  y <- .rho_np_y(rho_np = rho)
  A <- .calculate_A(y = y, n1 = n1, mu = mu, sigma = sigma, design = design)
  int <- hcubature(
    f = \(z_) {
      n2_ <- n2_extrapol(design, z_)
      x1_ <- .z_to_x(z_, n = n1, mu0 = mu0, sigma = sigma)
      (1 - pnorm((.rho_np_B(rho_np = rho, x1 = x1_, n1 = n1, n2 = n2_, mu0 = mu0, mu1 = mu1) - mu) * sqrt(n2_) / sigma)) *
        .f1(z1 = z_, n1 = n1, mu = mu, mu0 = mu0, sigma = sigma)
    },
    lowerLimit = c1f,
    upperLimit = c1e,
    vectorInterface = TRUE
  )$integral
  int + A
}
.rho_np_root <- function(gamma, x1, x2, mu0, mu1, sigma, design){
  uniroot(
    f = \(x) .p_np(x1 = x1, x2 = x2, mu = x, mu0 = mu0, mu1 = mu1, sigma = sigma, design = design) - gamma,
    interval = qnorm(c(.025, .975), mean = mu0, sd = sigma / sqrt(design@n1)),
    extendInt = "yes"
  )$root
}
#### (Median unbiased point estimator) (NP ordering).
.median_unbiased_np <- function(x1, x2, mu0, mu1, sigma, design){
  # (This should not really be used.)
  .rho_np_root(gamma = 0.5, x1 = x1, x2 = x2, mu0 = mu0, mu1 = mu1, sigma = sigma, design = design)
}
#### (Confidence interval) (NP ordering).
.confidence_interval_np <- function(alpha, x1, x2, mu0, mu1, sigma, design){
  # This should not really be used.
  c(.rho_np_root(gamma = alpha/2, x1 = x1, x2 = x2, mu0 = mu0, mu1 = mu1, sigma = sigma, design = design),
    .rho_np_root(gamma = 1-alpha/2, x1 = x1, x2 = x2, mu0 = mu0, mu1 = mu1, sigma = sigma, design = design))
}

### Stage-wise combination function ordering
#### P-value (SWCF ordering).
.p_swcf <- function(x1, x2, mu, mu0, sigma, design) {
  n1 <- design@n1
  c1f <- design@c1f
  c1e <- design@c1e
  z1 <- .x_to_z(x = x1, n = n1, mu0 = mu0, sigma = sigma)
  if ((z1 < c1f || z1 > c1e) && (x2 != -Inf)) {
    stop("z1 suggests early stopping but x2 is not -Inf.")
  }
  if ((z1 >= c1f  && z1 <= c1e) && (x2 == -Inf)) {
    stop("z1 suggests continuation but x2 is -Inf.")
  }
  ret <-
    if (z1 < c1f || z1 > c1e) {
      1 - pnorm( (x1 - mu)*sqrt(n1)/sigma )
    } else {
      n2 <- n2_extrapol(design, z1)
      cf_value <- (x2 - mu) * sqrt(n2) / sigma - c2_extrapol(design, (x1 - mu) * sqrt(n1) / sigma)
      A <- 1 - pnorm(c1e -mu*sqrt(n1)/sigma)
      int <- hcubature(
        f = \(z_) {
          n2_ <- n2_extrapol(design, z_)
          x1_ <- .z_to_x(z_, n = n1, mu0 = mu0, sigma = sigma)
          (1 - pnorm((.rho_swcf_B(cf_value = cf_value, x1 = x1_, n1 = n1, n2 = n2_, mu = mu) - mu) * sqrt(n2_) / sigma)) *
            .f1(z1 = z_, n1 = n1, mu = mu, mu0 = mu0, sigma = sigma)
        },
        lowerLimit = c1f,
        upperLimit = c1e,
        vectorInterface = TRUE
      )$integral
      A + int
    }
  ret
}
.rho_swcf_root <- function(gamma, x1, x2, mu0, sigma, design){
  uniroot(
    f = \(x) .p_swcf(x1 = x1, x2 = x2, mu = x, mu0 = mu0, sigma = sigma, design = design) - gamma,
    interval = qnorm(c(.025, .975), mean = mu0, sd = sigma / sqrt(design@n1)),
    extendInt = "yes"
  )$root
}
#### Median unbiased point estimator (SWCF ordering).
.median_unbiased_swcf <- function(x1, x2, mu0, sigma, design){
  .rho_swcf_root(gamma = 0.5, x1 = x1, x2 = x2, mu0 = mu0, sigma = sigma, design = design)
}
#### Confidence interval (SWCF ordering).
.confidence_interval_swcf <- function(alpha, x1, x2, mu0, sigma, design){
  c(.rho_swcf_root(gamma = alpha/2, x1 = x1, x2 = x2, mu0 = mu0, sigma = sigma, design = design),
    .rho_swcf_root(gamma = 1-alpha/2, x1 = x1, x2 = x2, mu0 = mu0, sigma = sigma, design = design))
}


### Naive CI
.naive_confidence_interval <- function(alpha, x1, x2, n1, n2, sigma){
  ret <-
    if (x2 == -Inf) {
      c(x1 - qnorm(1 - alpha/2)* sigma/sqrt(n1), x1 + qnorm(1 - alpha/2)* sigma/sqrt(n1))
    } else {
      x <- .x1_x2_to_x(x1 = x1, x2 = x2, n1 = n1, n2 = n2)
      c(x - qnorm(1 - alpha/2)* sigma/sqrt(n1 + n2), x + qnorm(1 - alpha/2)* sigma/sqrt(n1 + n2))
    }
  ret
}

### Repeated confidence interval.
.repeated_ci_l1 <- function(x1, sigma, design){
  x1 - design@c1e * sigma / sqrt(design@n1)
}
.repeated_ci_l2 <- function(x1, x2, mu0, sigma, design){
  n1 <- design@n1
  c1f <- design@c1f
  c1e <- design@c1e
  z1 <- .x_to_z(x = x1, n = n1, mu0 = mu0, sigma = sigma)
  n2 <- n2_extrapol(design, z1)
  lower_l <- x1 - c1e * sigma / sqrt(n1)
  upper_l <- x1 - c1f * sigma / sqrt(n1)
  lower_z2 <- (x2 - lower_l) * sqrt(n2) / sigma
  upper_z2 <- (x2 - upper_l) * sqrt(n2) / sigma
  minc2 <- c2_extrapol(design, c1e)
  maxc2 <- c2_extrapol(design, c1f)
  # Assumes c2 is monotonically decreasing
  if (maxc2 < lower_z2) {
    ret <- upper_l
  } else if (minc2 >upper_z2) {
    ret <- lower_l
  } else{
    ret <- uniroot(\(x) c2_extrapol(design, (x1 - x)*sqrt(n1)/sigma) - (x2 - x)*sqrt(n2)/sigma,
                   interval = c(lower_l, upper_l))$root
  }
  ret
}
.repeated_confidence_interval <- function(x1, x2, mu0, sigma, design) {
  if (x2 == -Inf){
    c(.repeated_ci_l1(x1 = x1, sigma = sigma, design = design),
      -.repeated_ci_l1(x1 = -x1, sigma = sigma, design = design))
  } else {
    c(.repeated_ci_l2(x1 = x1, x2 = x2, mu0 = mu0, sigma = sigma, design = design),
      -.repeated_ci_l2(x1 = -x1, x2 = -x2, mu0 = -mu0, sigma = sigma, design = design))
  }
}
