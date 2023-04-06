library(microbenchmark)

x1 <- runif(10, min = -10, max = 10)
h <- (designad@c1e - designad@c1f) / 2
x <- h * designad@x1_norm_pivots + (h + designad@c1f)
y <- designad@n2_pivots
ff1 <- stats::splinefun(x, y, method = "fmm")
ff2 <- stats::splinefun(x, y, method = "periodic")
ff3 <- stats::splinefun(x, y, method = "natural")
ff4 <- stats::splinefun(x, y, method = "monoH.FC")
ff5 <- stats::splinefun(x, y, method = "hyman")
fnew <- fastmonoH.FC(x, y)
microbenchmark(
  stats::splinefun(x, y, method = "fmm")(x1),
  ff1(x1),
  ff2(x1),
  ff3(x1),
  ff4(x1),
  ff5(x1),
  fnew(x1))

splinefun2 <- function (x, y = NULL, method = c("fmm", "periodic", "natural",
                                                "monoH.FC", "hyman"), ties = mean)
{
  x <- stats:::regularize.values(x, y, ties, missing(ties))
  y <- x$y
  x <- x$x
  nx <- length(x)
  if (is.na(nx))
    stop(gettextf("invalid value of %s", "length(x)"), domain = NA)
  if (nx == 0)
    stop("zero non-NA points")
  method <- match.arg(method)
  if (method == "periodic" && y[1L] != y[nx]) {
    warning("spline: first and last y values differ - using y[1L] for both")
    y[nx] <- y[1L]
  }
  if (method == "monoH.FC") {
    n1 <- nx - 1L
    dy <- y[-1L] - y[-nx]
    dx <- x[-1L] - x[-nx]
    Sx <- dy/dx
    m <- c(Sx[1L], (Sx[-1L] + Sx[-n1])/2, Sx[n1])
    m <- .Call(stats:::C_monoFC_m, m, Sx, PACKAGE = "stats")
    return(splinefunH02(x0 = x, y0 = y, m = m, dx = dx))
  }
  iMeth <- match(method, c("periodic", "natural", "fmm", "monoH.FC",
                           "hyman"))
  if (iMeth == 5L) {
    dy <- diff(y)
    if (!(all(dy >= 0) || all(dy <= 0)))
      stop("'y' must be increasing or decreasing")
  }
  z <- .Call(stats:::C_SplineCoef, min(3L, iMeth), x, y)
  return(z)
  if (iMeth == 5L)
    z <- spl_coef_conv(hyman_filter(z))
  rm(x, y, nx, method, iMeth, ties)
  function(x, deriv = 0L) {
    deriv <- as.integer(deriv)
    if (deriv < 0L || deriv > 3L)
      stop("'deriv' must be between 0 and 3")
    if (deriv > 0L) {
      z0 <- double(z$n)
      z[c("y", "b", "c")] <- switch(deriv, list(y = z$b,
                                                b = 2 * z$c, c = 3 * z$d), list(y = 2 * z$c,
                                                                                b = 6 * z$d, c = z0), list(y = 6 * z$d, b = z0,
                                                                                                           c = z0))
      z[["d"]] <- z0
    }
    res <- .splinefun(x, z)
    if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L]))
      res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
    res
  }
}
splinefunH02 <- function (x0, y0, m, dx = x0[-1L] - x0[-length(x0)])
{
  function(x, deriv = 0, extrapol = c("linear", "cubic")) {
    extrapol <- match.arg(extrapol)
    deriv <- as.integer(deriv)
    if (deriv < 0 || deriv > 3)
      stop("'deriv' must be between 0 and 3")
    i <- findInterval(x, x0, all.inside = (extrapol == "cubic"))
    if (deriv == 0)
      interp <- function(x, i) {
        h <- dx[i]
        t <- (x - x0[i])/h
        t1 <- t - 1
        h01 <- t * t * (3 - 2 * t)
        h00 <- 1 - h01
        tt1 <- t * t1
        h10 <- tt1 * t1
        h11 <- tt1 * t
        y0[i] * h00 + h * m[i] * h10 + y0[i + 1] * h01 +
          h * m[i + 1] * h11
      }
    else if (deriv == 1)
      interp <- function(x, i) {
        h <- dx[i]
        t <- (x - x0[i])/h
        t1 <- t - 1
        h01 <- -6 * t * t1
        h10 <- (3 * t - 1) * t1
        h11 <- (3 * t - 2) * t
        (y0[i + 1] - y0[i])/h * h01 + m[i] * h10 + m[i +
                                                       1] * h11
      }
    else if (deriv == 2)
      interp <- function(x, i) {
        h <- dx[i]
        t <- (x - x0[i])/h
        h01 <- 6 * (1 - 2 * t)
        h10 <- 2 * (3 * t - 2)
        h11 <- 2 * (3 * t - 1)
        ((y0[i + 1] - y0[i])/h * h01 + m[i] * h10 +
            m[i + 1] * h11)/h
      }
    else interp <- function(x, i) {
      h <- dx[i]
      h01 <- -12
      h10 <- 6
      h11 <- 6
      ((y0[i + 1] - y0[i])/h * h01 + m[i] * h10 + m[i +
                                                      1] * h11)/h
    }
    if (extrapol == "linear" && any(iXtra <- (iL <- (i ==
                                                     0)) | (iR <- (i == (n <- length(x0)))))) {
      r <- x
      if (any(iL))
        r[iL] <- if (deriv == 0)
          y0[1L] + m[1L] * (x[iL] - x0[1L])
      else if (deriv == 1)
        m[1L]
      else 0
      if (any(iR))
        r[iR] <- if (deriv == 0)
          y0[n] + m[n] * (x[iR] - x0[n])
      else if (deriv == 1)
        m[n]
      else 0
      ini <- !iXtra
      r[ini] <- interp(x[ini], i[ini])
      r
    }
    else {
      interp(x, i)
    }
  }
}

ff4 <- splinefun2(x, y, method = "monoH.FC")
fastmonoH.FC <- function(x, y=NULL, ties = mean) {
  x <- stats:::regularize.values(x, y, ties, missing(ties))
  y <- x$y
  x <- x$x
  nx <- length(x)
  if (is.na(nx))
    stop(gettextf("invalid value of %s", "length(x)"), domain = NA)
  if (nx == 0)
    stop("zero non-NA points")
  n1 <- nx - 1L
  dy <- y[-1L] - y[-nx]
  dx <- x[-1L] - x[-nx]
  Sx <- dy/dx
  m <- c(Sx[1L], (Sx[-1L] + Sx[-n1])/2, Sx[n1])
  m <- .Call(stats:::C_monoFC_m, m, Sx, PACKAGE = "stats")
  p0 <- y[-length(y)]
  p1 <- y[-1L]
  ddx <-  x[-1L] - x[-length(x)]
  m0 <- m[-length(m)] * ddx
  m1 <- m[-1L] * ddx
  b <- m0 / ddx
  c <- (3* (p1 - p0) - 2*m0 -m1) / ddx^2
  d <- (2 * (p0 - p1) + m0 + m1) / ddx^3
  z <- list(
    method=3L,
    n=n1,
    x=x,
    y=y,
    b=b,
    c=c,
    d=d
  )
  rm(x, y, nx, ties, n1, dy, dx, Sx, m, p0, p1, ddx, m0, m1, b, c, d)
  function(x, deriv = 0L) {
    deriv <- as.integer(deriv)
    if (deriv < 0L || deriv > 3L)
      stop("'deriv' must be between 0 and 3")
    if (deriv > 0L) {
      z0 <- double(z$n)
      z[c("y", "b", "c")] <- switch(deriv, list(y = z$b,
                                                b = 2 * z$c, c = 3 * z$d), list(y = 2 * z$c,
                                                                                b = 6 * z$d, c = z0), list(y = 6 * z$d, b = z0,
                                                                                                           c = z0))
      z[["d"]] <- z0
    }
    res <- stats:::.splinefun(x, z)
    if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L]))
      res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
    res
  }
}


myspline3 <- function(x) {
  i <- findInterval(x, x0)
  dx <- x - x0[i]
  y0[i] + dx * (b[i] + dx * (c[i] + dx * d[i]))
}
myspline2 <- function(x) {
  i <- findInterval(x, x0)
  dx = x0[-1L] - x0[-length(x0)]
  h <- dx[i]
  t <- (x - x0[i])/h
  t1 <- t - 1
  h01 <- t * t * (3 - 2 * t)
  h00 <- 1 - h01
  tt1 <- t * t1
  h10 <- tt1 * t1
  h11 <- tt1 * t
  y0[i] * h00 + h * m[i] * h10 + y0[i + 1] * h01 +
    h * m[i + 1] * h11
}
myspline <- fastmonoH.FC(x,y)
myspline(x1)
ff4(x1, extrapol = "c")


myspline2(a)
other_spline <- splinefun2(x, y, method = "fmm")

x0 <- get("x0", envir = environment(ff4))
y0 <- get("y0", envir = environment(ff4))

p <- y0
m <- get("m", envir = environment(ff4))
p0 <- p[-length(p)]
p1 <- p[-1L]
ddx <-  x0[-1L] - x0[-length(x0)]
m0 <- m[-length(m)] * ddx
m1 <- m[-1L] * ddx
b <- m0 / ddx
c <- (3* (p1 - p0) - 2*m0 -m1) / ddx^2
d <- (2 * (p0 - p1) + m0 + m1) / ddx^3

a <- runif(10, min=x0[1], max=x0[length(x0)])
y + x(b[i] + x (c + x (d)))

microbenchmark(
 stats::splinefun(x, y, method = "fmm")(x1),
 ff(x1, extrapol = "cubic"),
 ff(x1, extrapol = "linear"),
 ff1(x1),
 ff2(x1),
 ff3(x1),
 ff4(x1),
 ff5(x1),
 myspline(x1))




find_d2_1(designad, 0.2, 0.2, 1, FALSE, 0.5)
















