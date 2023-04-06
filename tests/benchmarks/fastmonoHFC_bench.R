library(microbenchmark)
x10 <- runif(10)
y10 <- runif(10)
x100 <- runif(100)
y100 <- runif(100)
x1000 <- runif(1000)
y1000 <- runif(1000)

x10s <- sort(x10)
y10s <- sort(y10)

microbenchmark(
  oldspline <- splinefun(x10, y10, method = "monoH.FC"),
  newspline <- fastmonoH.FC(x10, y10)
)
microbenchmark(
  oldspline <- splinefun(x100, y100, method = "monoH.FC"),
  newspline <- fastmonoH.FC(x10, y10)
)
microbenchmark(
  oldspline <- splinefun(x100, y100, method = "monoH.FC"),
  newspline <- fastmonoH.FC(x100, y100)
)
microbenchmark(
  oldspline <- splinefun(x1000, y1000, method = "monoH.FC"),
  newspline <- fastmonoH.FC(x1000, y1000)
)

microbenchmark(
  oldspline <- splinefun(x10s, y10s, method = "monoH.FC"),
  newspline <- fastmonoH.FC(x10s, y10s)
)




a10 <- rnorm(10)
a100 <- rnorm(100)
a1000 <- rnorm(1000)



oldspline <- splinefun(x10, y10, method = "monoH.FC")
newspline <- fastmonoH.FC(x10, y10, extrapol = "linear")
oldspline(a10, extrapol = "linear")
newspline(a10)

splinecoef <- get("z", envir = environment(newspline))


microbenchmark(
 oldspline(a10),
 newspline(a10)
)

microbenchmark(
  oldspline(a100),
  newspline(a100)
)

microbenchmark(
  oldspline(a1000),
  newspline(a1000)
)




