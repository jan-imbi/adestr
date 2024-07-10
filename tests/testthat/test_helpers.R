test_that("Transformation convenience functions work",
{
  expect_equal(t_to_z(z_to_t(1.234, 1.3, 2.1, 23), 1.3, 2.1, 23), 1.234)
  expect_equal(smean_to_t(1.2, 3.4, 45, FALSE), 1.2/sqrt(3.4) * sqrt(45))
  expect_equal(smeans_to_smean(1, 2, 30, 30), 1.5)
  expect_equal(svar_to_se(2.3, 40, FALSE), sqrt(2.3/40))
  expect_equal(t_test(1.2, 3.4, 50, FALSE), 1.2/sqrt(3.4)*sqrt(50))
  expect_equal(pool_svar1_svar2(1.2, 3.4, 50, 60, FALSE), (49 * 1.2 + 59 * 3.4) / (108))
})



test_that("Generation of example statistics doesn't produce an error",
{
  expect_error(get_example_statistics(), NA)
  expect_error(get_statistics_from_paper (), NA)
})
