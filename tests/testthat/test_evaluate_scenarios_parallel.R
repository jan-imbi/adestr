test_that("evaluate_scenarios_parallel works for one scenario without error",
          {
            expect_error(
              evaluate_scenarios_parallel(
                score_lists = list(c(MSE(), OverestimationProbability())),
                estimator_lists =  list(c(SampleMean(), FirstStageSampleMean())),
                data_distribution_lists = list(c(Normal(FALSE), Normal(TRUE))),
                design_lists =  list(c(designad, designgs)),
                mu_lists = list(c(-1, 0, 1)),
                sigma_lists = list(1)
              )
              ,
              NA
            )
          })
test_that("evaluate_scenarios_parallel works for two scenarios without error",
          {
            expect_error(
              evaluate_scenarios_parallel(
                score_lists = list(c(MSE(), OverestimationProbability()), c(Width())),
                estimator_lists =  list(c(SampleMean(), FirstStageSampleMean()), c(NaiveCI())),
                data_distribution_lists = list(c(Normal(FALSE), Normal(TRUE)), c(Normal(FALSE))),
                design_lists =  list(c(designad, designgs), c(designgs)),
                mu_lists = list(c(-1, 0, 1), c(0)),
                sigma_lists = list(1, 1)
              )
              ,
              NA
            )
          })
test_that("evaluate_scenarios_parallel works for two scenarios without error",
          {
            expect_error(
              evaluate_scenarios_parallel(
                score_lists = MSE(),
                estimator_lists =  c(SampleMean(), FirstStageSampleMean()),
                data_distribution_lists = Normal(),
                design_lists =  list(c(designad, designgs)),
                mu_lists = c(-1, 0, 1),
                sigma_lists = 1
              )
              ,
              NA
            )
          })
test_that("conversion of non-lists to lists works properly",
          {
            expect_error(
              evaluate_scenarios_parallel(
                score_lists = MSE(),
                estimator_lists =  c(SampleMean(), FirstStageSampleMean()),
                data_distribution_lists = Normal(),
                design_lists =  list(c(designad, designgs)),
                mu_lists = c(-1, 0, 1),
                sigma_lists = 1
              ),
              NA
            )
          })


