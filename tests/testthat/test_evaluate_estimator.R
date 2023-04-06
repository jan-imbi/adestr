test_that("MSE of sample mean can be calculated without error.",
          {
            expect_error(
              evaluate_estimator(
                score = MSE(),
                estimator = SampleMean(),
                data_distribution = Normal(),
                design = designad,
                mu = c(0.3),
                sigma = 1,
                exact = FALSE
              )
              ,
              NA
            )
          })








