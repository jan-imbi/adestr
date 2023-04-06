datadist <- Normal(two_armed = FALSE)
H_0 <- PointMassPrior(.0, 1)
H_1 <- PointMassPrior(.4, 1)
ess <- ExpectedSampleSize(datadist, H_1)
power <- Power(datadist, H_1)
toer  <- Power(datadist, H_0)
ad_initialD <- get_initial_design(
  theta = .4,
  alpha = .025,
  beta  = .2,
  type_design  = "two-stage",
  dist  = Normal(two_armed = FALSE)
)
opt <- minimize(
  ess,
  subject_to(
    power >= 0.8,
    toer  <= .025
  ),
  ad_initialD
)
design <- opt$design
# nVecs <- c(64L, 128L, 256L, 512L, 1024L)
nVecs <- c(100L)
# methods <- c("hcubature", "pcubature", "cuhre", "divonne", "suave", "vegas")
# divonne too slow??
methods <- c("hcubature", "pcubature", "cuhre", "suave", "vegas")
library(cubature)


nested_integrals <- function(smean2, n1, n2, mu, sigma, two_armed){
  sapply(smean2,
         \(x) integrate(
           f1_1,
           -Inf,
           Inf,
           rel.tol = 1e-5,
           smean2 = x,
           n1 = n1,
           n2 = n2,
           mu = mu,
           sigma = sigma,
           two_armed = two_armed
         )$value)
}
start_time <- Sys.time()
invisible(replicate(1e1, integrate(nested_integrals, -Inf, Inf, rel.tol = 1e-5, n1 = 30, n2 = 30, mu = 0, sigma = 1, two_armed = FALSE)))
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
invisible(replicate(1e1, hcubature(
  mf1_1,
  lowerLimit = c(-Inf, -Inf),
  upperLimit = c(Inf, Inf),
  tol = 1e-7,
  vectorInterface = TRUE,
  n1 = 30, n2 = 30, mu = 0, sigma = 1, two_armed = FALSE
  )))
end_time <- Sys.time()
end_time - start_time






results <- data.frame()
for (nVec in nVecs){
  for (method in methods) {
    t_start <- Sys.time()
    res <- int_design_2(design, method = method, nVec = nVec,  mu=0, sigma = 1, two_armed = FALSE)
    t_end <- Sys.time()
    results <- rbind(results, data.frame(
      method = method,
      nVec = nVec,
      elapsed_time = t_end-t_start,
      integral = res$futility_integral$integral + res$continuation_integral$integral + res$efficacy_integral$integral
    ))
  }
}
results


start_time <- Sys.time()
invisible(replicate(1e4, integrate(dnorm, qnorm(1e-12, lower.tail = TRUE), qnorm(1e-12, lower.tail = FALSE), subdivisions = 1e3) ))
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
invisible(replicate(1e4, hcubature(dnorm, qnorm(1e-12, lower.tail = TRUE), qnorm(1e-12, lower.tail = FALSE), tol = .Machine$double.eps^0.25, vectorInterface = TRUE ) ))
end_time <- Sys.time()
end_time - start_time




relTol <- 1e-05
absTol <- 1e-12
maxEval <- 10^6
two_armed <- FALSE
mu <- 0
sigma <- 1
n1 <- n1(design, round=FALSE)
se1 <- (1L + two_armed) *  sigma/sqrt(n1)
df1 <- (1L + two_armed) * (n1 - 1L)
n2_min <- n2(design, design@c1f, round=FALSE)
n2_max <- n2(design, design@c1e, round=FALSE)
n_min <- n1 + n2_min
n_max <- n1 + n2_max
df_min <- (1L + two_armed) * (n_min - 1L)
df_max <- (1L + two_armed) * (n_max - 1L)
se_min <- (1L + two_armed) * sigma/sqrt(n_min)
se_max <- (1L + two_armed) * sigma/sqrt(n_max)
mf1 <- \(x, n1, mu, sigma, two_armed) {
  n <- n1 + n2(design, x[3,], round=FALSE)
  mtf1_2(x, n, n1, mu, sigma, two_armed)
  print(mtf1_2(x, n, n1, mu, sigma, two_armed))
}
min_norm <- min(qnorm(absTol, mean = mu, sd=se_min, lower.tail = TRUE),
                qnorm(absTol, mean = mu, sd=se_max, lower.tail = TRUE))
max_norm <- min(qnorm(absTol, mean = mu, sd=se_min, lower.tail = FALSE),
                qnorm(absTol, mean = mu, sd=se_max, lower.tail = FALSE))
min_chisq <- min(qchisq(absTol, df = df_min, lower.tail = TRUE) * sigma^(2L) / df_min,
                 qchisq(absTol, df = df_max, lower.tail = TRUE) * sigma^(2L) / df_max)
max_chisq <- max(qchisq(absTol, df = df_min, lower.tail = FALSE) * sigma^(2L) / df_min,
                 qchisq(absTol, df = df_max, lower.tail = FALSE) * sigma^(2L) / df_max)

library(cubature)

start_t <- Sys.time()
a <- cubintegrate(
  f = mf1,
  lower = c(min_norm,
            min_chisq,
            design@c1f,
            qchisq(absTol, df = df1, lower.tail = TRUE) * sigma^(2L) / df1),
  upper = c(max_norm,
            max_chisq,
            design@c1e,
            qchisq(absTol, df = df1, lower.tail = FALSE) * sigma^(2L) / df1),
  fDim = 1L,
  method = "hcubature",
  nVec = 1073741824L,
  relTol = 1e-4,
  absTol = absTol,
  maxEval =10^9,
  n1 = n1, mu = mu, sigma = sigma, two_armed = two_armed
)
end_t <- Sys.time()
end_t - start_t






library(cubature)
# cubintegrate(\(x) 1, 0, 1, method="vegas",  stateFile = "a.txt")


cuhre(\(x) {1}, lowerLimit = 0, upperLimit = 1,
      flags = list(verbose = 0L, final = 1L, smooth = 0L, keep_state = 1L, load_state = 0L,
                                                    level = 0L),
      stateFile = "a.txt")




cubintegrate(
  f = mf1,
  lower = c(min_norm,
            min_chisq,
            design@c1f,
            qchisq(absTol, df = df1, lower.tail = TRUE) * sigma^(2L) / df1),
  upper = c(max_norm,
            max_chisq,
            design@c1e,
            qchisq(absTol, df = df1, lower.tail = FALSE) * sigma^(2L) / df1),
  fDim = 1L,
  method = "vegas",
  relTol = 1e-2,
  absTol = 1e-2,
  maxEval = maxEval,
  nVec = 100L,
  flags = list(keep_state = TRUE),
  stateFile = "a.txt",
  n1, mu, sigma, two_armed
)


















