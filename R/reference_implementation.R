.x1_x2_to_x <- function(x1, x2, n1, n2){
  (n1 * x1 + n2 * x2) / (n1 + n2)
}
.n_sigma_to_se <- function(n, sigma) {
  sigma/sqrt(n)
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
.f1 <- function(z1, n1, mu, mu0, sigma) {
  dnorm(z1 + (mu0 - mu) * sqrt(n1)/sigma)
}
.f2 <- function(z1, z2, n1, n2, mu, mu0, sigma){
  .f1(z1 = z1, n1 = n1, mu = mu, mu0 = mu0, sigma = sigma) * dnorm(z2 + (mu0 - mu) * sqrt(n2)/sigma)
}

.rho_ml1 <- function(x1){
  x1
}
.rho_ml2 <- function(x1, x2, n1, n2){
  .x1_x2_to_x(x1 = x1, x2 = x2, n1 = n1, n2 = n2)
}
.rho_lr1 <- function(x1, n1, mu) {
  (x1 - mu) * sqrt(n1)
}
.rho_lr2 <- function(x1, x2, n1, n2, mu) {
  (.x1_x2_to_x(x1 = x1, x2 = x2, n1 = n1, n2 = n2) - mu) * sqrt(n1 + n2)
}
.rho_st1 <- function(x1, n1, mu) {
  n1 * x1 - n1 * mu
}
.rho_st2 <- function(x1, x2, n1, n2, mu) {
  n1 * x1 + n2 * x2 - (n1 + n2) * mu
}

.rho_np1 <- function(x1, n1, mu0, mu1) {
  (2*x1*(mu0 - mu1) - mu0 - mu1)*n1
}
.rho_np2 <- function(x1, x2, n1, n2, mu0, mu1) {
  (2*.x1_x2_to_x(x1 = x1, x2 = x2, n1 = n1, n2 = n2)*(mu0 - mu1) - mu0 - mu1)*(n1 + n2)
}
.rho_swcf1 <- function(x1, n1, mu0, sigma){
  pnorm(.x_to_z(x = x1, n = n1, mu0 = mu0, simga = sigma))
}
.rho_swcf2 <- function(x1, x2, n1, n2, mu, sigma, design){
  1 + .swcf(x1 = x1, x2 = x2, n1 = n1, n2 = n2, mu = mu, sigma = sigma, design = design)
}
.rho_swcf3 <- function(x1, n1, mu0, sigma){
  2 + pnorm(.x_to_z(x = x1, n = n1, mu0 = mu0, simga = sigma))
}
.swcf <- function(x1, x2, n1, n2, mu, sigma, design){
  pnorm(.x_to_z(x = x2, n = n2, mu0 = mu, sigma = sigma) - c2_extrapol(design, .x_to_z(x = x1, n = n1, mu0 = mu, sigma = sigma)))
}

# This doesn't work :(
.mle_cdf <- function(y, mu, mu0, sigma, design) {
  n1 <- design@n1
  yz <- .x_to_z(x = y, n = n1, mu0 = mu0, sigma = sigma)
  domain <-
  if (yz < design@c1f)
    c(-Inf, yz)
  else if (yz > design@c1e)
    c(design@c1e, yz)
  else
    c(-Inf, design@c1f)
  early_stopping_int <- hcubature(
    .f1,
    vectorInterface = TRUE,
    lowerLimit = domain[1],
    upperLimit = domain[2],
    n1 = n1,
    mu = mu,
    mu0 = mu0,
    sigma = sigma
    )$integral
  continuation_int <- hcubature(
    \(x, n1, n2, mu, mu0, sigma) .f2(z1 = x[1,,drop=FALSE],
                                     z2 = x[2,,drop=FALSE],
                                     n1 = n1,
                                     n2 = n2,
                                     mu = mu,
                                     mu0 = mu0,
                                     sigma = sigma
                                     ),
    vectorInterface = TRUE,
    lowerLimit = c(design@c1f, -Inf),
    upperLimit = c(design@c1e, (y-mu0)*(n1 + n2) - z1*sqrt(n1)/sigma),
    n1 = n1,
    mu = mu,
    mu0 = mu0,
    sigma = sigma
  )$integral
  continuation_int
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
    n1 = n1,
    mu = mu,
    mu0 = mu0,
    sigma = sigma,
    design = design
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
    n1 = n1, mu = mu, mu0 = mu0, sigma = sigma)$integral

  efficacy <- hcubature(
    f = \(x, n1, mu, mu0, sigma){
      .z_to_x(z = x[1,,drop=FALSE], n = n1, mu0 = mu0, sigma = sigma) *
        .f1(z1 = x[1,,drop=FALSE], n1 = n1, mu = mu, mu0 = mu0, sigma = sigma)
    },
    lowerLimit = c(design@c1e),
    upperLimit = c(Inf),
    vectorInterface = TRUE,
    n1 = n1,
    mu = mu,
    mu0 = mu0,
    sigma = sigma
  )$integral
  continuation <- hcubature(
    f = \(x, n1, mu, mu0, sigma, design){
      n2 <- n2_extrapol(design, x[1,,drop=FALSE])
      .z1_z2_to_x(
        z1 = x[1, , drop = FALSE],
        z2 = x[2, , drop = FALSE],
        n1 = n1,
        n2 = n2,
        mu0 = mu0,
        sigma = sigma
      ) *
        .f2(
          z1 = x[1, , drop = FALSE],
          z2 = x[2, , drop = FALSE],
          n1 = n1,
          n2 = n2,
          mu = mu,
          mu0 = mu0,
          sigma = sigma
        )
    },
    lowerLimit = c(design@c1f,-Inf),
    upperLimit = c(design@c1e, Inf),
    vectorInterface = TRUE,
    n1 = n1,
    mu = mu,
    mu0 = mu0,
    sigma = sigma,
    design = design
  )$integral
  efficacy + futility + continuation
}
.mle_expectation_derivative <- function(mu, mu0, sigma, design){
  n1 <- design@n1
  futility <- hcubature(
    f = \(x, n1, mu, mu0, sigma){
      .z_to_x(z = x[1,,drop=FALSE], n = n1, mu0 = mu0, sigma = sigma) *
        (x[1,,drop=FALSE] + (mu0 - mu) * sqrt(n1)/sigma) * sqrt(n1)/sigma *
        .f1(z1 = x[1,,drop=FALSE], n1 = n1, mu = mu, mu0 = mu0, sigma = sigma)
    },
    lowerLimit = c(-Inf),
    upperLimit = c(design@c1f),
    vectorInterface = TRUE,
    n1 = n1,
    mu = mu,
    mu0 = mu0,
    sigma = sigma
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
    n1 = n1,
    mu = mu,
    mu0 = mu0,
    sigma = sigma
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
    n1 = n1,
    mu = mu,
    mu0 = mu0,
    sigma = sigma,
    design = design
  )$integral
  efficacy + futility + continuation
}

.bias_reduced <- function(mle, mu, mu0, sigma, iterations, design){
  estimate <- mle
  browser()
  for (i in seq_len(iterations)) {
    expectation_mle <- .mle_expectation(mu = estimate, mu0 = mu0, sigma = sigma, design = design)
    derivative_mle <- .mle_expectation_derivative(mu = estimate, mu0 = mu0, sigma = sigma, design = design)
    estimate <- estimate +  ((mle - estimate) - (expectation_mle - mle)) / (1 + derivative_mle)
  }
  estimate
}


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










