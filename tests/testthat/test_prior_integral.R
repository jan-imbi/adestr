test_that("Calculating MSE wrt. a non-degenerate-prior is roughly the same as wrt. to a point-prior.",
{
  expect_equal(
    evaluate_estimator(
      score = MSE(),
      estimator = SampleMean(),
      data_distribution = Normal(),
      use_full_twoarm_sampling_distribution = FALSE,
      design = get_example_design(),
      mu = NormalPrior(mu=0, sigma=1),
      sigma = 1
    )@results$MSE,
    evaluate_estimator(
      score = MSE(),
      estimator = SampleMean(),
      data_distribution = Normal(),
      use_full_twoarm_sampling_distribution = FALSE,
      design = get_example_design(),
      mu = 0,
      sigma = 1
    )@results$MSE,
    tolerance=1e-1
    )

  expect_equal(
    evaluate_estimator(
      score = MSE(),
      estimator = SampleMean(),
      data_distribution = Normal(),
      use_full_twoarm_sampling_distribution = FALSE,
      design = get_example_design(),
      mu = UniformPrior(min = -.1, max = .1),
      sigma = 1
    )@results$MSE,
    evaluate_estimator(
      score = MSE(),
      estimator = SampleMean(),
      data_distribution = Normal(),
      use_full_twoarm_sampling_distribution = FALSE,
      design = get_example_design(),
      mu = 0,
      sigma = 1
    )@results$MSE,
    tolerance=1e-2
  )
}
)

test_that("get_pdf works.",
{
  expect_equal(
    get_pdf(NormalPrior())(0), dnorm(0)
    )
  expect_equal(
    get_pdf(UniformPrior())(0), dunif(0, min=-1, max=1)
  )
})




