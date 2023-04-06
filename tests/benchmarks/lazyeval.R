library(microbenchmark)
f1 <- function(x) return(x)
f2 <- function(x, ...) return(x)
f3 <- function(x, ...) {
  args <- list(...)
  x
}

a <- runif(100000)

microbenchmark(
  f1(3),
  f1(x = 3),
  f2(x = 3),
  f2(x = 3, y=a),
  f2(x = 3, y=mean(a), z=cos(a), yy=mean(cos(a)), zz = mean(sin(a))),
  f3(x = 3, y=mean(a), z=cos(a), yy=mean(cos(a)), zz = mean(sin(a)))
)




