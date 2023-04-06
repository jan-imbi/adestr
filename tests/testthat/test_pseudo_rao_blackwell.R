test_that("Pseudo Rao-Blackwell estimator has less variance than Rao-Blackwell estimator",
          {
            prb <- evaluate_estimator(score = Variance(),
                                      estimator = PseudoRaoBlackwell(),
                                      design = designad,
                                      data_distribution = Normal(FALSE),
                                      mu = 0.3,
                                      sigma = 1
            )
            rb <- evaluate_estimator(score = Variance(),
                                     estimator = RaoBlackwell(),
                                     design = designad,
                                     data_distribution = Normal(FALSE),
                                     mu = 0.3,
                                     sigma = 1
            )
            expect_lt(
              prb@results$Variance,
              rb@results$Variance)
            })
