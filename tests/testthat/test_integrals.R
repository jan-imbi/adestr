test_that("integral over sample space is equal to 1 for case: known variance, one-armed",
          {
            expect_equal(
              int_kv(
                design = designad,
                mu = .3,
                sigma = 2.1,
                two_armed = FALSE
              )$overall_integral$integral,
              1,
              tolerance=1e-5
            )
          })
test_that("integral over sample space is equal to 1 for case: known variance, one-armed (exact=TRUE)",
          {
            expect_equal(
              int_kv(
                design = designad,
                mu = .3,
                sigma = 2.1,
                two_armed = FALSE,
                exact = TRUE
              )$overall_integral$integral,
              1,
              tolerance=1e-5
            )
          })
test_that("integral over sample space is equal to 1 for case: unknown variance, one-armed",
          {
            expect_equal(
              int_uv(
                design = designad,
                mu = 0.3,
                sigma = 2.1,
                two_armed = FALSE)$overall_integral$integral,
              1,
              tolerance=1e-3
            )
          })
test_that("integral over sample space is equal to 1 for case: unknown variance, one-armed (exact=TRUE)",
          {
            expect_equal(
              int_uv(
                design = designad,
                mu = 0.3,
                sigma = 2.1,
                two_armed = FALSE,
                exact=TRUE)$overall_integral$integral,
              1,
              tolerance=1e-3
            )
          })
test_that("integral over sample space is equal to 1 for case: known variance, two-armed",
          {
            expect_equal(
              int_kv(
                design = designad,
                mu = .3,
                sigma = 2.1,
                two_armed = TRUE)$overall_integral$integral,
              1,
              tolerance=1e-5
            )
          })
test_that("integral over sample space is equal to 1 for case: known variance, two-armed (exact=TRUE)",
          {
            expect_equal(
              int_kv(
                design = designad,
                mu = .3,
                sigma = 2.1,
                two_armed = TRUE,
                exact=TRUE)$overall_integral$integral,
              1,
              tolerance=1e-5
            )
          })
test_that("integral over sample space is equal to 1 for case: unknown variance, two-armed",
          {
            expect_equal(
              int_uv(
                design = designad,
                mu = .3,
                sigma = 2.1,
                two_armed = TRUE)$overall_integral$integral,
              1,
              tolerance=1e-5
            )
          })
test_that("integral over sample space is equal to 1 for case: unknown variance, two-armed (exact=TRUE)",
          {
            expect_equal(
              int_uv(
                design = designad,
                mu = .3,
                sigma = 2.1,
                two_armed = TRUE,
                exact=TRUE)$overall_integral$integral,
              1,
              tolerance=1e-5
            )
          })
test_that("integral over sample space is equal to 1 for case: known variance, two-armed, full sampling distribution",
          {
            expect_equal(
              int_kv_full(
                design = designad,
                mu = .3,
                sigma = 2.1)$overall_integral$integral,
              1,
              tolerance=1e-5
            )
          })
test_that("integral over sample space is equal to 1 for case: known variance, two-armed, full sampling distribution (exact=TRUE)",
          {
            expect_equal(
              int_kv_full(
                design = designad,
                mu = .3,
                sigma = 2.1,
                exact=TRUE)$overall_integral$integral,
              1,
              tolerance=1e-5
            )
          })
test_that("integral over sample space is equal to 1 for case: unknown variance, two-armed, full sampling distribution",
          {
            expect_equal(
              int_uv_full(
                design = designad,
                mu = .3,
                sigma = 2.1)$overall_integral$integral,
              1,
              tolerance=3e-2
            )
          })
test_that("integral over sample space is equal to 1 for case: unknown variance, two-armed, full sampling distribution, exact n2",
          {
            expect_equal(
              int_uv_full(
                design = designad,
                mu = .3,
                sigma = 2.1,
                exact = TRUE)$overall_integral$integral,
              1,
              tolerance=5e-2
            )
          })

