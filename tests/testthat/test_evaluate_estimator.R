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


test_that("evaluate_estimator works on lists without error.",
          {
            expect_error(
              evaluate_estimator(
                score = list(Bias(), Variance()),
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



test_that("Bias of sample mean can be calculated without error.",
          {
            expect_error(
              evaluate_estimator(
                score = Bias(),
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

test_that("Coverage of NaiveCI can be calculated without error.",
          {
            expect_error(
              evaluate_estimator(
                score = Coverage(),
                estimator = NaiveCI(),
                data_distribution = Normal(),
                design = designad,
                mu = c(0.3),
                sigma = 1,
                exact = FALSE
              )
              ,
              NA
            )
            expect_error(
              evaluate_estimator(
                score = SoftCoverage(),
                estimator = NaiveCI(),
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



test_that("TestAgreement of NaiveCI can be calculated",
          {
            expect_lt(
              evaluate_estimator(
                score = TestAgreement(),
                estimator = NaiveCI(),
                data_distribution = Normal(),
                design = designad,
                mu = c(0.3),
                sigma = 1,
                exact = FALSE
              )@results$`Agreement with test decision`
              ,
              1
            )
          })

test_that("TestAgreement of NaivePValue can be calculated",
          {
            skip_on_cran()
            expect_lt(
              evaluate_estimator(
                score = TestAgreement(),
                estimator = NaivePValue(),
                data_distribution = Normal(),
                design = designad,
                mu = c(0.3),
                sigma = 1,
                exact = FALSE
              )@results$`Agreement with test decision`
              ,
              1
            )
          })


test_that("Centrality of SampleMean wrt. NaiveCI can be calculated",
          {
            expect_equal(
              evaluate_estimator(
                score = Centrality(interval = NaiveCI()),
                estimator = SampleMean(),
                data_distribution = Normal(),
                design = designad,
                mu = c(0.3),
                sigma = 1,
                exact = FALSE
              )@results$Centrality
              ,
              0
            )
          })








