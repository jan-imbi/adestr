#' @importFrom grDevices xy.coords
Myregularize.values <- function (x, y, ties, warn.collapsing = TRUE, na.rm = TRUE)
{
  x <- xy.coords(x, y, setLab = FALSE)
  y <- x$y
  x <- x$x
  keptNA <- FALSE
  nx <- if (any(na <- is.na(x) | is.na(y))) {
    ok <- !na
    if (na.rm) {
      x <- x[ok]
      y <- y[ok]
      length(x)
    }
    else {
      keptNA <- TRUE
      sum(ok)
    }
  }
  else {
    length(x)
  }
  if (!identical(ties, "ordered")) {
    ordered <- if (is.function(ties) || is.character(ties))
      FALSE
    else if (is.list(ties) && length(T <- ties) == 2L &&
             is.function(T[[2]])) {
      ties <- T[[2]]
      identical(T[[1]], "ordered")
    }
    else stop("'ties' is not \"ordered\", a function, or list(<string>, <function>)")
    if (!ordered && is.unsorted(if (keptNA)
      x[ok]
      else x)) {
      o <- order(x)
      x <- x[o]
      y <- y[o]
    }
    if (length(ux <- unique(x)) < nx) {
      if (warn.collapsing)
        warning("collapsing to unique 'x' values")
      y <- as.vector(tapply(y, match(x, x), ties))
      x <- ux
      stopifnot(length(y) == length(x))
      if (keptNA)
        ok <- !is.na(x)
    }
  }
  list(x = x, y = y, keptNA = keptNA, notNA = if (keptNA) ok)
}
fastmonoH.FC_function <- function(x, y=NULL, ties = mean, extrapol = c("linear", "cubic")) {
  x <- Myregularize.values(x, y, ties, missing(ties))
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
  m <- .Call(MymonoFC_m, m, Sx, PACKAGE = "adestr")
  p0 <- y[-length(y)]
  p1 <- y[-1L]
  ddx <-  x[-1L] - x[-length(x)]
  m0 <- m[-length(m)] * ddx
  m1 <- m[-1L] * ddx
  b <- m0 / ddx
  c <- (3* (p1 - p0) - 2*m0 -m1) / ddx^2
  d <- (2 * (p0 - p1) + m0 + m1) / ddx^3
  extrapol <- match.arg(extrapol)
  iextrapol <- match(extrapol, c("linear", "cubic"))
  z <- list(
    method=iextrapol,
    n=n1,
    x=x,
    y=y,
    b=b,
    c=c,
    d=d,
    ml=m[1L],
    mr=m[length(m)]
  )
  rm(x, y, nx, ties, n1, dy, dx, Sx, m, p0, p1, ddx, m0, m1, b, c, d, extrapol, iextrapol)
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
    res <- .Call(MySplineEval, x, z, PACKAGE = "adestr")
    if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L]))
      res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
    res
  }
}

fastmonoH.FC_coefficients <- function(x, y=NULL, ties = mean, extrapol = c("linear", "cubic")){
  x <- Myregularize.values(x, y, ties, missing(ties))
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
  m <- .Call(MymonoFC_m, m, Sx, PACKAGE = "adestr")
  p0 <- y[-length(y)]
  p1 <- y[-1L]
  ddx <-  x[-1L] - x[-length(x)]
  m0 <- m[-length(m)] * ddx
  m1 <- m[-1L] * ddx
  b <- m0 / ddx
  c <- (3* (p1 - p0) - 2*m0 -m1) / ddx^2
  d <- (2 * (p0 - p1) + m0 + m1) / ddx^3
  extrapol <- match.arg(extrapol)
  iextrapol <- match(extrapol, c("linear", "cubic"))
  z <- list(
    method=iextrapol,
    n=n1,
    x=x,
    y=y,
    b=b,
    c=c,
    d=d,
    ml=m[1L],
    mr=m[length(m)]
  )
  rm(x, y, nx, ties, n1, dy, dx, Sx, m, p0, p1, ddx, m0, m1, b, c, d, extrapol, iextrapol)
  z
}
fastmonoH.FC_evaluate <- function(x, z, deriv = 0L) {
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
  res <- .Call(MySplineEval, x, z, PACKAGE = "adestr")
  if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L]))
    res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
  res
}
