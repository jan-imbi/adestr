library(cubature)
design <- get_example_design()

# View design characteristics
plot(design)

# Density integrates to 1
continuation <-
  hcubature(
    f = \(z) .f(z[1], z[2], mu = 0, mu0 = 0, sigma = 1, design = design),
    lowerLimit = c(design@c1f, -Inf),
    upperLimit = c(design@c1e, Inf),
    vectorInterface = FALSE
  )
futility <-
  hcubature(
    f = \(z) .f(z[1], -Inf, mu = 0, mu0 = 0, sigma = 1, design = design),
    lowerLimit = c(-Inf),
    upperLimit = c(design@c1f),
    vectorInterface = FALSE
  )
efficacy <-
  hcubature(
    f = \(z) .f(z[1], -Inf, mu = 0, mu0 = 0, sigma = 1, design = design),
    lowerLimit = c(design@c1e),
    upperLimit = c(Inf),
    vectorInterface = FALSE
  )
continuation$integral + futility$integral + efficacy$integral

# Plot distribution of the MLE
x <- seq(-0.5, 1.5, 0.01)
y <- sapply(x, .mle_pdf, mu = 0.2, mu0 = 0, sigma = 1, design = design)
plot(x, y, type = "l")

# Plot bias of the MLE
mu <- seq(-0.5, 1, 0.05)
mle_bias <- function(mu, mu0, sigma, design) {
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
  efficacy + futility + continuation - mu
}
y <- sapply(mu, mle_bias, mu0 = 0, sigma = 1, design = design)
plot(mu, y, type = "l")

# Example: Coverage of the LR ordering confidence interval
lr_ordering_ci_coverage <- function(alpha, mu, mu0, sigma, design) {
  n1 <- design@n1
  futility <- hcubature(
    f = \(z, n1, mu, mu0, sigma){
      x1 <- .z_to_x(z = z, n = n1, mu0 = mu0, sigma = sigma)
      ci <- sapply(x1, .confidence_interval_lr, x2 = -Inf, mu0 = mu0, sigma = sigma, design = design, alpha = alpha)
      (ci[1,,drop=FALSE] < mu  & mu < ci[2,,drop=FALSE]) *
        .f1(z1 = z[1,,drop=FALSE], n1 = n1, mu = mu, mu0 = mu0, sigma = sigma)
    },
    lowerLimit = c(-Inf),
    upperLimit = c(design@c1f),
    vectorInterface = TRUE,
    n1 = n1, mu = mu, mu0 = mu0, sigma = sigma
  )$integral
  efficacy <- hcubature(
    f = \(z, n1, mu, mu0, sigma){
      x1 <- .z_to_x(z = z, n = n1, mu0 = mu0, sigma = sigma)
      ci <- sapply(x1, .confidence_interval_lr, x2 = -Inf, mu0 = mu0, sigma = sigma, design = design, alpha = alpha)
      (ci[1,,drop=FALSE] < mu  & mu < ci[2,,drop=FALSE]) *
        .f1(z1 = z[1,,drop=FALSE], n1 = n1, mu = mu, mu0 = mu0, sigma = sigma)
    },
    lowerLimit = c(design@c1e),
    upperLimit = c(Inf),
    vectorInterface = TRUE,
    n1 = n1, mu = mu, mu0 = mu0, sigma = sigma
  )$integral
  continuation <- hcubature(
    f = \(z, n1, mu, mu0, sigma, design){
      n2 <- n2_extrapol(design, z[1,,drop=FALSE])
      x1 <- .z_to_x(z = z[1,,drop=FALSE], n = n1, mu0 = mu0, sigma = sigma)
      x2 <- .z_to_x(z = z[2,,drop=FALSE], n = n2, mu0 = mu0, sigma = sigma)
      ci <- mapply(x1 = x1, x2 = x2,
                   FUN = .confidence_interval_lr, MoreArgs = list(mu0 = mu0, sigma = sigma,  design = design, alpha = alpha))
      (ci[1,,drop=FALSE] < mu  & mu < ci[2,,drop=FALSE]) *
        .f2(z1 = z[1, , drop = FALSE], z2 = z[2, , drop = FALSE], n1 = n1, n2 = n2, mu = mu, mu0 = mu0, sigma = sigma)
    },
    lowerLimit = c(design@c1f,-Inf),
    upperLimit = c(design@c1e, Inf),
    vectorInterface = TRUE,
    tol = 1e-3, # default 1e-5 would take very long.
    n1 = n1, mu = mu, mu0 = mu0, sigma = sigma, design = design
  )$integral
  efficacy + futility + continuation
}
# This may take a minute
lr_ordering_ci_coverage(alpha = 0.05, mu = 0, mu0 = 0, sigma = 1, design = design)




