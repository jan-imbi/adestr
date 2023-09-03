datadist <- Normal(two_armed = FALSE)
H_0 <- PointMassPrior(.0, 1)
# The power is calculated under the point hypothesis mu=0.4.
H_1 <- PointMassPrior(.4, 1)
ess_H0 <- ExpectedSampleSize(datadist, H_0)
ess_H1 <- ExpectedSampleSize(datadist, H_1)
toer  <- Power(datadist, H_0)
power <- Power(datadist, H_1)
# Here, the intial design parameters from which the optimization will start are set.
initial_ad <- get_initial_design(
  theta = .4,
  alpha = .025,
  beta  = .2,
  type_design  = "two-stage",
  dist  = datadist,
  type_n2 = "linear_decreasing",
  cf = 0,
  ce = 2.1
)
evaluate(toer, initial_ad)
initial_ad@tunable[["c1e"]] <- FALSE
# initial_ad@tunable[["c1f"]] <- FALSE
# Here, the parameters for the adaptive design are optimized.
designad <- minimize(
  ess_H1,
  subject_to(
    power >= 0.8,
    toer  <= .025
  ),
  initial_ad
)$design
plot(designad)
designad <- cache_design_splines(designad)
evaluate(toer, designad)

adestr:::plot_p(LikelihoodRatioOrderingPValue(), Normal(two_armed = FALSE),
                designad, 0, 1, boundary_color = scales::hue_pal()(5)[3], subdivisions = 200L)


adestr:::plot_p(StagewiseCombinationFunctionOrderingPValue(), Normal(two_armed = FALSE),
                designad, 0, 1, boundary_color = scales::hue_pal()(5)[3], subdivisions = 200L)


