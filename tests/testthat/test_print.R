test_that("Console output for various classes doesn't produce errors.",
{
  set.seed(321)
  dat <- data.frame(
    endpoint = rnorm(sum(c(56, 56, 47, 47)), mean = rep(c(.3, 0, .3, 0), c(56, 56, 47, 47))),
    group = factor(rep(c("ctl", "trt", "ctl", "trt"), c(56,56,47,47))),
    stage = rep(c(1L, 2L), c(56*2, 47*2))
  )
  expect_error(
    capture.output(show(designad)),
    NA)
  expect_error(
    capture.output(show(evaluate_estimator(
      score = MSE(),
      estimator = SampleMean(),
      data_distribution = Normal(),
      design = designad,
      mu = c(0.3),
      sigma = 1,
      exact = FALSE
    ))),
    NA)
  expect_error(
    capture.output(show(analyze(
      data = dat,
      statistics = list(FirstStageSampleMean(), SampleMean(), NaiveCI()),
      data_distribution = Normal(TRUE),
      sigma = 1,
      design = get_example_design(TRUE)
    ))),
    NA)

}
)

