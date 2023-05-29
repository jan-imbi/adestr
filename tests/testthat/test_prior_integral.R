evaluate_estimator(
  score = MSE(),
  estimator = SampleMean(),
  data_distribution = Normal(),
  use_full_twoarm_sampling_distribution = FALSE,
  design = get_example_design(),
  mu = NormalPrior(mu=0, sigma=1),
  sigma = 1
  )

