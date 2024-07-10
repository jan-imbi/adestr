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





