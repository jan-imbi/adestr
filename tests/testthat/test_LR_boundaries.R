# test_that("implicit rejection boundary from Neyman-Pearson test matches c2",
#           {
#             p2 <- get_stagewise_estimators(NeymanPearsonOrderingPValue(0, 0.4),
#                                      Normal(FALSE),
#                                      FALSE,
#                                      designad,
#                                      1,
#                                      FALSE)[[2]]
#             expect_equal(
#               implied_c2(
#                 designad,
#                 adoptr:::scaled_integration_pivots(designad),
#                 p2,
#                 1,
#                 FALSE,
#                 0.025
#               ),
#               designad@c2_pivots,
#               tolerance = 1e-3
#             )
#           })
