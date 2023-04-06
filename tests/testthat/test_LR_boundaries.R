test_that("implicit rejection boundary from Neyman-Pearson test matches c2",
          {
            expect_equal(
              c2_np(
                designad,
                adoptr:::scaled_integration_pivots(designad),
                0,
                0,
                0.4,
                1,
                FALSE,
                0.025
              ),
              designad@c2_pivots,
              tolerance = 1e-3
            )
          })
