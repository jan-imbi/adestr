library(microbenchmark)
library(mvtnorm)

x <- matrix(c(runif(100),runif(100)),nrow=2)
xx <- t(x)

.INTEGRAL_CONST <- 1/(2 * pi *  sqrt(1-0.5^2))
.EXP_CONST <- -1/(2*(1-0.5^2))
Sigma <- matrix(c(1, 0.5, 0.5, 1),ncol=2)
# mdbinorm <- function(x, mu, sigma) {
#   xT <- (x[1L,,drop=FALSE] - mu)
#   .INTEGRAL_CONST / sigma^2 * exp(.EXP_CONST/sigma^2*(xT * (xT-x[2L,,drop=FALSE]) + x[2L,,drop=FALSE]^2))
# }
mdbinorm2 <- function(x, mu, sigma) .INTEGRAL_CONST / sigma^2 * exp(.EXP_CONST/sigma^2*((x[1L,,drop=FALSE]-mu) * ((x[1L,]-mu)-x[2L,,drop=FALSE]) + x[2L,,drop=FALSE]^2))
mdbinorm3 <- function(x, mu, sigma) {
  xT <- (x[1L,] - mu)
  xC <- x[2L,]
  .INTEGRAL_CONST / sigma^2 * exp(.EXP_CONST/sigma^2*(xT * (xT-xC) + xC^2))
}

microbenchmark(
  mdbinorm(x, 0, 1),
  mdbinorm2(x, 0, 1),
  mdbinorm3(x, 0, 1),
  dmvnorm(xx, c(0,0), Sigma, checkSymmetry = FALSE, log = FALSE)
)

mldbinorm <- function(x, muT, muC, sigma) {
  xT <- (x[1L,] - muT)
  yT <- (x[2L,] - muT)
  log(.INTEGRAL_CONST / sigma^2) + .EXP_CONST/sigma^2*(xT * (xT-yT) + yT^2)
}

a <- rnorm(100)
b <- rnorm(100)

x <- matrix(c(a,b),nrow=2)
xx <- t(x)

sigma <- 2

microbenchmark(
  mdbinorm(x, 0, 0, 1),
  mldbinorm(x, 0, 0, 1)
)


sqrt(1-.5^2)
microbenchmark(
dmvnorm(xx, c(0,0), Sigma, checkSymmetry = FALSE, log = TRUE),
1/(2*pi*sigma^2*sqrt(1-.5^2))*exp(-1/(2*(1-.5^2)) * ( (x[1L,]/sigma)^2 -2*.5*(x[1L,]/sigma)*(x[2L,]/sigma) + (x[2L,]/sigma)^2 )),
dnorm(x[1L,]) * dnorm(x[2L,])*exp(-(x[1L,]/sigma) * (x[2L,]/sigma)),
mdbinorm(x, 0, 0, 1),
mldbinorm(x, 0, 0, 1)
)

y <- .3
ddx <- function(y) {
  delta <- .0001
  y_1 <- hcubature(
    mdbinorm,
    c(qnorm(1e-7), qnorm(1e-7)),
    c(y-delta/2, 0),
    vectorInterface = TRUE,
    mu = 0,
    sigma = 1
  )$integral
  y_2 <-  hcubature(
    mdbinorm,
    c(qnorm(1e-7), qnorm(1e-7)),
    c(y+delta/2, 0),
    vectorInterface = TRUE,
    mu = 0,
    sigma = 1
  )$integral
  (y_2 - y_1)/delta
}

hcubature(\(x, mu, sigma) mdbinorm(rbind(x, y), mu, sigma),
          c(qnorm(1e-7)),
          c(0),
          vectorInterface = TRUE,
          mu = 0,
          sigma = 1)

prob_boundary_inf <- if (absError==0) min(.Machine$double.eps^.25, tol) else min(.Machine$double.eps^.25, absError)



hcubature(\(x){
  n2 <- n2_extrapol(design, x[1L,]/se1)
  mdbinorm(rbind(x[1L,], (y*(n1 + n2) - n2 * x[3L,])/n1 ), mu, sigma) *
    mf0_1((y*(n1 + n2) - n1 * x[2L,])/n2, mu, sigma, FALSE) *
    mf0_1(x[4L,], mu, sigma, FALSE) *
    (n1 + n2)^2/ (n1*n2)
},
c(design@c1f * se1,
  qnorm(prob_boundary_inf, mean = mu, sd = se1_onearm, lower.tail = TRUE),
  qnorm(prob_boundary_inf, mean = mu, sd = se2_onearm, lower.tail = TRUE),
  qnorm(prob_boundary_inf, mean = mu, sd = se2_onearm, lower.tail = TRUE)),
c(design@c1e * se1,
  qnorm(prob_boundary_inf, mean = mu, sd = se1_onearm, lower.tail = FALSE),
  qnorm(prob_boundary_inf, mean = mu, sd = se2_onearm, lower.tail = FALSE),
  qnorm(prob_boundary_inf, mean = mu, sd = se2_onearm, lower.tail = FALSE)),
vectorInterface = TRUE,
mu = 0,
sigma = 1)












