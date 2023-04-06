test_that(
  "MLMSE score can be used to optimize designs, and the resuling desing has a higher first-stage sample size.",
  {
    H_0 <- PointMassPrior(.0, 1)
    H_1 <- PointMassPrior(.4, 1)
    ess <- ExpectedSampleSize(Normal(FALSE), H_1)
    power <- Power(Normal(FALSE), H_1)
    toer  <- Power(Normal(FALSE), H_0)
    initial_design <- get_initial_design(
      theta = .4,
      alpha = .025,
      beta  = .2,
      type_design  = "two-stage",
      dist  = Normal(FALSE),
      order = 7L
    )
    mlmse <- MLMSE(mu = c(0, .3, .6))
    cs <- composite({
      ess + 100 * mlmse
    })
    opt <- minimize(cs,
                    subject_to(power >= 0.8,
                               toer  <= .025),
                    initial_design)
    expect_gt(opt$design@n1,
              designad@n1)
  }
)
