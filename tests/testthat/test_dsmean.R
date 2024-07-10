test_that("density of MLE sums up to one (normal distribution, one-armed)",
          {
            expect_equal(
              dsmean(
                Normal(two_armed = FALSE),
                designad,
                .x <- seq(-2, 2, .h <- .01),
                0,
                1,
                exact = FALSE,
                combine_components = TRUE
              ) |> sum() * .h,
              1,
              tolerance=1e-2
            )
          })
test_that("density of MLE sums up to one (normal distribution, one-armed) (exact=TRUE)",
          {
            expect_equal(
              dsmean(
                Normal(two_armed = FALSE),
                designad,
                .x <- seq(-2, 2, .h <- .01),
                0,
                1,
                exact = TRUE,
                combine_components = TRUE
              ) |> sum() * .h,
              1,
              tolerance=1e-2
            )
          })

test_that("density of MLE sums up to one (t distribution, one-armed)",
          {
            skip_on_cran()
            expect_equal(
              dsmean(
                Student(two_armed = FALSE),
                designad,
                .x <- seq(-2, 2, .h <- .001),
                0,
                1,
                exact = FALSE,
                combine_components = TRUE
              ) |> sum() * .h,
              1,
              tolerance = 1e-2
            )
          })

test_that("density of MLE sums up to one (t distribution, one-armed) (exact=TRUE)",
          {
            skip_on_cran()
            expect_equal(
              dsmean(
                Student(two_armed = FALSE),
                designad,
                .x <- seq(-2, 2, .h <- .001),
                0,
                1,
                exact = TRUE,
                combine_components = TRUE
              ) |> sum() * .h,
              1,
              tolerance = 1e-2
            )
          })


test_that("density of MLE sums up to one (normal distribution, two-armed, treatment group)",
          {
            skip_on_cran()
            expect_equal(dsmeanT(Normal(),
                                 designad,
                                 .x <- seq(-2, 2, .h <- .01),
                                 0,
                                 1,
                                 exact = FALSE) |> sum() * .h,
                         1,
                         tolerance = 1e-2)
          })

test_that("density of MLE sums up to one (normal distribution, two-armed, treatment group) (exact=TRUE)",
          {
            skip_on_cran()
            expect_equal(dsmeanT(Normal(),
                                 designad,
                                 .x <- seq(-2, 2, .h <- .01),
                                 0,
                                 1,
                                 exact = TRUE) |> sum() * .h,
                         1,
                         tolerance = 1e-2)
          })

test_that("density of MLE sums up to one (t distribution, two-armed, treatment group)",
          {
            skip_on_cran()
            expect_equal(dsmeanT(Student(),
                                 designad,
                                 .x <- seq(-2, 2, .h <- .01),
                                 0,
                                 1,
                                 exact = FALSE) |> sum()*.h,
                         1,
                         tolerance = 1e-2)
          })

test_that("density of MLE sums up to one (t distribution, two-armed, treatment group) (exact=TRUE)",
          {
            skip_on_cran()
            expect_equal(dsmeanT(Student(),
                                 designad,
                                 .x <- seq(-2, 2, .h <- .1),
                                 0,
                                 1,
                                 exact = TRUE) |> sum()*.h,
                         1,
                         tolerance = 1e-2)
          })




