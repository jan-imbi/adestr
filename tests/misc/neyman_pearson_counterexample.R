

datadist <- Normal(two_armed = FALSE)
H_0 <- PointMassPrior(.0, 1)
H_1 <- PointMassPrior(.4, 1)
H_01 <- PointMassPrior(c(0, .4), c(.5, .5))
ess <- ExpectedSampleSize(datadist, H_0)
toer  <- Power(datadist, H_0)
power <- Power(datadist, H_1)
cp <- ConditionalPower(datadist, H_1)
initial_ad <- get_initial_design(
  theta = .4,
  alpha = .025,
  beta  = .2,
  type_design  = "two-stage",
  dist  = datadist,
  type_n2 = "linear_increasing"
)
optad <- minimize(
  ess,
  subject_to(
    power >= 0.8,
    toer  <= .025
  ),
  initial_ad
  # opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-06, maxeval = 40000)
)
plot(optad$design)

design <- adestr:::cache_design_splines(optad$design)


stagewise_estimators <- get_stagewise_estimators(estimator = NeymanPearsonOrderingPValue(0, 0.4),
                                                 data_distribution =  datadist,
                                                 use_full_twoarm_sampling_distribution = FALSE,
                                                 design = design, sigma = 1, exact = FALSE)
p1 <- stagewise_estimators[[1L]]
p2 <- stagewise_estimators[[2L]]



implied_c2(design, seq(design@c1f, design@c1e, .1), p2, 1, FALSE, 0.025) -
c2(design, seq(design@c1f, design@c1e, .1))

