devtools::load_all(".")

# Calculate all the designs
single_stage <- rpact::getSampleSizeMeans(groups=1, normalApproximation = TRUE, alternative = 0.4)
datadist <- Normal(two_armed = FALSE)
H_0 <- PointMassPrior(.0, 1)
H_1 <- PointMassPrior(.4, 1)
H_01 <- PointMassPrior(c(0, .4), c(.5, .5))
ess_H0 <- ExpectedSampleSize(datadist, H_0)
ess_H1 <- ExpectedSampleSize(datadist, H_1)
ess_H01 <- ExpectedSampleSize(datadist, H_01)
power <- Power(datadist, H_1)
toer  <- Power(datadist, H_0)
cp <- ConditionalPower(datadist, H_1)
initial_gs <- get_initial_design(
  theta = .4,
  alpha = .025,
  beta  = .2,
  type_design  = "group-sequential",
  dist  = datadist
)
initial_gs@tunable[c("c1f", "n1")] <- c(FALSE)
initial_gs@n1 <- 10
initial_ad <- get_initial_design(
  theta = .4,
  alpha = .025,
  beta  = .2,
  type_design  = "two-stage",
  dist  = datadist
)
initial_ad@tunable[c("c1f", "n1")] <- c(FALSE)
initial_ad@n1 <- 10
optad <- minimize(
  ess_H1,
  subject_to(
    power >= 0.8,
    toer  <= .025
  ),
  initial_ad,
  opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-05, maxeval = 40000)
)
designad_noattr <- optad$design
designad <- cache_design_splines(designad_noattr)

t0 <- Sys.time()
evaluate_estimator(ProbOverestimation(),
                   MedianUnbiased(),
                   Normal(two_armed=FALSE),
                   designad_noattr,
                   mu=0.2,
                   sigma=1,
                   exact=FALSE)
t1 <- Sys.time()
print(sprintf("Time for non-cached version: %s", format(t1 - t0)))


h <- (designad@c1e - designad@c1f) / 2
attr(designad, "c2_cache") <-
  stats::splinefun(h * designad@x1_norm_pivots + (h + designad@c1f),
                   designad@c2_pivots, method = "monoH.FC")
attr(designad, "n2_cache") <-
  stats::splinefun(h * designad@x1_norm_pivots + (h + designad@c1f),
                   designad@n2_pivots, method = "monoH.FC")
t0 <- Sys.time()
evaluate_estimator(ProbOverestimation(),
                   BiasReduced(),
                   Normal(two_armed=FALSE),
                   designad,
                   mu=0.2,
                   sigma=1,
                   exact=FALSE)
t1 <- Sys.time()
print(sprintf("Time for cached version: %s", format(t1 - t0)))

attr(designad, "c2_cache") <-
  fastmonoH.FC(h * designad@x1_norm_pivots + (h + designad@c1f),
                   designad@c2_pivots)
attr(designad, "n2_cache") <-
  fastmonoH.FC(h * designad@x1_norm_pivots + (h + designad@c1f),
               designad@n2_pivots)
t0 <- Sys.time()
evaluate_estimator(ProbOverestimation(),
                   BiasReduced(),
                   Normal(two_armed=FALSE),
                   designad,
                   mu=0.2,
                   sigma=1,
                   exact=FALSE)
t1 <- Sys.time()
print(sprintf("Time for cached version with fast splines: %s", format(t1 - t0)))










