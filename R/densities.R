# Constants
.LSQRT2PIINV <- log(1/sqrt(2*pi))
.TLSQRT2PIINV <- 2*.LSQRT2PIINV
.LTPISQRTHALFINV <- log(1/(2*pi*sqrt(.5)))
.TLTPISQRTHALFINV <- 2*.LTPISQRTHALFINV
.SQRT2 <- sqrt(2)

# Known variance
logf1_kv <- function(z1, n1, mu, sigma, two_armed) .LSQRT2PIINV - 0.5* (z1 - mu / sigma * sqrt(n1 / (1L + two_armed)) )^2
logf2_kv <- function(z1, z2, n1, n2, mu, sigma, two_armed) {
  musig <- mu/sigma
  .TLSQRT2PIINV - 0.5 * ((z1 - musig *  sqrt(n1 / (1L + two_armed)))^2 + (z2 - musig * sqrt(n2 / (1L + two_armed)) )^2)
}
logf1_kv_full <- function(z1D, z1C, n1, mu, sigma) {
  t1D_shift <- z1D -  mu * sqrt(n1/2L) / sigma
  .LTPISQRTHALFINV - t1D_shift * (t1D_shift + .SQRT2*z1C) - z1C^2
}
logf2_kv_full <- function(z1D, z1C, t2D, t2C, n1, n2, mu, sigma) {
  musig <- mu/(sigma * sqrt(2L))
  t1D_shift <- z1D -  musig * sqrt(n1)
  t2D_shift <- t2D -  musig * sqrt(n2)
  .TLTPISQRTHALFINV - t1D_shift * (t1D_shift + .SQRT2*z1C) - z1C^2 - t2D_shift * (t2D_shift + .SQRT2*t2C) - t2C^2
}
f1_kv <- function(z1, n1, mu, sigma, two_armed) exp(logf1_kv(z1, n1, mu, sigma, two_armed))
f2_kv <- function(z1, z2, n1, n2, mu, sigma, two_armed) exp(logf2_kv(z1, z2, n1, n2, mu, sigma, two_armed))
f1_kv_full <- function(z1D, z1C, n1, mu, sigma) exp(logf1_kv_full(z1D, z1C, n1, mu, sigma))
f2_kv_full <- function(z1D, z1C, t2D, t2C, n1, n2, mu, sigma) exp(logf2_kv_full(z1D, z1C, t2D, t2C, n1, n2, mu, sigma))
mf1_kv <- function(x, n1, mu, sigma, two_armed) f1_kv(x[1L,,drop=FALSE], n1, mu, sigma, two_armed)
mf2_kv <- function(x, n1, n2, mu, sigma, two_armed) f2_kv(x[1L,,drop=FALSE], x[2L,,drop=FALSE], n1, n2, mu, sigma, two_armed)
mf1_kv_full <- function(x, n1, mu, sigma) f1_kv_full(x[1L,,drop=FALSE], x[2L,,drop=FALSE], n1, mu, sigma)
mf2_kv_full <- function(x, n1, n2, mu, sigma) f2_kv_full(x[1L,,drop=FALSE], x[2L,,drop=FALSE], x[3L,,drop=FALSE], x[4L,,drop=FALSE], n1, n2, mu, sigma)
mlogf1_kv <- function(x, n1, mu, sigma, two_armed) logf1_kv(x[1L,,drop=FALSE], n1, mu, sigma, two_armed)
mlogf2_kv <- function(x, n1, n2, mu, sigma, two_armed) logf2_kv(x[1L,,drop=FALSE], x[2L,,drop=FALSE], n1, n2, mu, sigma, two_armed)
mlogf1_kv_full <- function(x, n1, mu, sigma) logf1_kv_full(x[1L,,drop=FALSE], x[2L,,drop=FALSE], n1, mu, sigma)
mlogf2_kv_full <- function(x, n1, n2, mu, sigma) logf2_kv_full(x[1L,,drop=FALSE], x[2L,,drop=FALSE], x[3L,,drop=FALSE], x[4L,,drop=FALSE], n1, n2, mu, sigma)

# Unknown variance
logf1_uv <- function(t1, svar1, n1, mu, sigma, two_armed) {
  v <- sigma^2
  varrat <- svar1/v
  sdrat <- sqrt(varrat)
  df <- (1L + two_armed) * (n1 - 1L)
  .LSQRT2PIINV - 0.5* (t1*sdrat - mu / sigma * sqrt(n1 / (1L + two_armed)))^2 + dchisq(varrat * df, df, log=TRUE) + log(df/v * sdrat)
}
logf2_uv <- function(t1, svar1, t2, svar2, n1, n2, mu, sigma, two_armed) {
  v <- sigma^2
  musig <- mu/sigma
  varrat1 <- svar1/v
  sdrat1 <- sqrt(varrat1)
  df1 <- (1L + two_armed) * (n1 - 1L)
  varrat2 <- svar2/v
  sdrat2 <- sqrt(varrat2)
  df2 <- (1L + two_armed) * (n2 - 1L)
  .TLSQRT2PIINV -
    0.5* ((t1*sdrat1 - musig * sqrt(n1 / (1L + two_armed)))^2 + (t2*sdrat2 - musig * sqrt(n2 / (1L + two_armed)))^2) +
    dchisq(varrat1 * df1, df1, log=TRUE) +
    dchisq(varrat2 * df2, df2, log=TRUE) +
    log(df1 * df2/(v^2)  * sdrat1 * sdrat2)
}
logf1_uv_full <- function(t1D, t1C, svar1, n1, mu, sigma) {
  musig <- mu/(.SQRT2 * sigma)
  v <- sigma^2
  varrat <- svar1/v
  sdrat <- sqrt(varrat)
  df <- 2L * (n1 - 1L)
  t1D_shift <- t1D*sdrat -  musig * sqrt(n1)
  t1C_mult <- t1C *sdrat
  .LTPISQRTHALFINV - t1D_shift * (t1D_shift + .SQRT2*t1C_mult) - t1C_mult^2 + dchisq(varrat * df, df, log=TRUE) + log(df/v * sdrat^2)
}
logf2_uv_full <- function(t1D, t1C, svar1, t2D, t2C, svar2, n1, n2, mu, sigma) {
  v <- sigma^2
  musig <- mu/(.SQRT2 * sigma)
  varrat1 <- svar1/v
  sdrat1 <- sqrt(varrat1)
  df1 <- 2L * (n1 - 1L)
  varrat2 <- svar2/v
  sdrat2 <- sqrt(varrat2)
  df2 <- 2L * (n2 - 1L)
  t1D_shift <- t1D*sdrat1 -  musig * sqrt(n1)
  t1C_mult <- t1C *sdrat1
  t2D_shift <- t2D*sdrat2 -  musig * sqrt(n2)
  t2C_mult <- t2C *sdrat2
  .TLTPISQRTHALFINV -
    t1D_shift * (t1D_shift + .SQRT2*t1C_mult) - t1C_mult^2 - t2D_shift * (t2D_shift + .SQRT2*t2C_mult) - t2C_mult^2 +
    dchisq(varrat1 * df1, df1, log=TRUE) +
    dchisq(varrat2 * df2, df2, log=TRUE) +
    log(df1 * df2/(v^2)  * sdrat1^2 * sdrat2^2)
}
f1_uv <- function(t1, svar1, n1, mu, sigma, two_armed) exp(logf1_uv(t1, svar1, n1, mu, sigma, two_armed))
f2_uv <- function(t1, svar1, t2, svar2, n1, n2, mu, sigma, two_armed) exp(logf2_uv(t1, svar1, t2, svar2, n1, n2, mu, sigma, two_armed))
f1_uv_full <- function(t1D, t1C, svar1, n1, mu, sigma) exp(logf1_uv_full(t1D, t1C, svar1, n1, mu, sigma))
f2_uv_full <- function(t1D, t1C, svar1, t2D, t2C, svar2, n1, n2, mu, sigma) exp(logf2_uv_full(t1D, t1C, svar1, t2D, t2C, svar2, n1, n2, mu, sigma))
mf1_uv <- function(x, n1, mu, sigma, two_armed) f1_uv(x[1L,,drop=FALSE], x[2L,,drop=FALSE], n1, mu, sigma, two_armed)
mf2_uv <- function(x, n1, n2, mu, sigma, two_armed) f2_uv(x[1L,,drop=FALSE], x[2L,,drop=FALSE], x[3L,,drop=FALSE], x[4L,,drop=FALSE], n1, n2, mu, sigma, two_armed)
mf1_uv_full <- function(x, n1, mu, sigma) f1_uv_full(x[1L,,drop=FALSE], x[2L,,drop=FALSE], x[3L,,drop=FALSE], n1, mu, sigma)
mf2_uv_full <- function(x, n1, n2, mu, sigma) f2_uv_full(x[1L,,drop=FALSE], x[2L,,drop=FALSE], x[3L,,drop=FALSE], x[4L,,drop=FALSE], x[5L,,drop=FALSE], x[6L,,drop=FALSE],n1, n2, mu, sigma)
mlogf1_uv <- function(x, n1, mu, sigma, two_armed) logf1_uv(x[1L,,drop=FALSE], x[2L,,drop=FALSE], n1, mu, sigma, two_armed)
mlogf2_uv <- function(x, n1, n2, mu, sigma, two_armed) logf2_uv(x[1L,,drop=FALSE], x[2L,,drop=FALSE], x[3L,,drop=FALSE], x[4L,,drop=FALSE], n1, n2, mu, sigma, two_armed)
mlogf1_uv_full <- function(x, n1, mu, sigma) logf1_uv_full(x[1L,,drop=FALSE], x[2L,,drop=FALSE], x[3L,,drop=FALSE], n1, mu, sigma)
mlogf2_uv_full <- function(x, n1, n2, mu, sigma) logf2_uv_full(x[1L,,drop=FALSE], x[2L,,drop=FALSE], x[3L,,drop=FALSE], x[4L,,drop=FALSE], x[5L,,drop=FALSE], x[6L,,drop=FALSE],n1, n2, mu, sigma)


# Densities with respect to sample mean (legacy code)
logf1_kv_smean <- function(smean1, n1, mu, sigma) {
  se1inv <- sqrt(n1)/sigma
  .TLSQRT2PIINV + log(se1inv) - 0.5* ((smean1 - mu) * se1inv)^2
}
logf2_smean <- function(smean1, smean2, n1, n2, mu, sigma) {
  se1inv <- sqrt(n1)/sigma
  se2inv <- sqrt(n2)/sigma
  .TLSQRT2PIINV + log(se1inv*se2inv) - 0.5* (((smean1 - mu) * se1inv)^2 + ((smean2 - mu) * se2inv)^2)
}


