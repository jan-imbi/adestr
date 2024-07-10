test_that("densities agree with reference implementation.",
{
  expect_equal(.f1(z1 = 1.5, n1 = 30, mu = 0, mu0 = 0, sigma = 1),
               f1_kv(z1 = 1.5, n1 = 30, mu = 0, sigma = 1, two_armed = FALSE))
  expect_equal(.f2(z1 = 1.5, z2 = 1, n1 = 30, n2 = 45, mu = 0, mu0 = 0, sigma = 1),
               f2_kv(z1 = 1.5, z2 = 1, n1 = 30, n2 = 45, mu = 0, sigma = 1, two_armed = FALSE))
  expect_equal(.f(z1 = 1.5, z2 = 1, mu = 0, mu0 = 0, sigma = 1, designad),
               f2_kv(z1 = 1.5, z2 = 1, n1 = 30, n2 = 45, mu = 0, sigma = 1, two_armed = FALSE))
})
test_that("MLE is correctly implemented.",
{
  expect_equal(.mle(1,2, 10, 20), (1*10 + 2*20)/30)
  expect_equal(.fixed_weighted_mle(1,2, 1), 1)

})
test_that("MLE densities agree with reference implementation.",
          {
            expect_equal(.mle_pdf(0.5, 0, 0, 1, designad),
                         dsmean(
                           Normal(two_armed = FALSE),
                           designad,
                           .5,
                           0,
                           1,
                           exact = FALSE,
                           combine_components = TRUE
                         ), tolerance = 1e-4)
            expect_equal(.mle_pdf(1.5, 1, 1, 1, designad),
                         dsmean(
                           Normal(two_armed = FALSE),
                           designad,
                           .5,
                           0,
                           1,
                           exact = FALSE,
                           combine_components = TRUE
                         ), tolerance = 1e-4)
          })
test_that("mle CDF works.",
{
  expect_equal(.mle_cdf(5, 0, 0, 1, designad), 1)
  expect_equal(.mle_cdf(-5, 0, 0, 1, designad), 0)
}
)


test_that("pseudo Rao-Blackwell estimator agrees with reference implementation.",
          {
            prb <- get_stagewise_estimators(PseudoRaoBlackwell(), Normal(FALSE), FALSE, designad, 1, FALSE)
            expect_equal(.pseudo_rao_blackwell(.4, .4, .0, 1, designad),
                         prb$g2(designad, .4, .4, designad@n1, n2_extrapol(designad, smean_to_z(.4, designad@n1, 1, FALSE)),
                                0, 1, FALSE),
                         tolerance = 1e-3)
          })

test_that("Rao-Blackwell estimator agrees with reference implementation.",
          {
            rb <- get_stagewise_estimators(RaoBlackwell(), Normal(FALSE), FALSE, designad, 1, FALSE)
            expect_equal(.rao_blackwell(.4, .4, .0, 1, designad),
                         rb$g2(.4, .4, designad@n1, n2_extrapol(designad, smean_to_z(.4, designad@n1, 1, FALSE)),
                                0, 1, FALSE),
                         tolerance = 1e-3)
          })
test_that("AWSM (min peak var) estimator agrees with reference implementation.",
          {
            awsm <- get_stagewise_estimators(MinimizePeakVariance(), Normal(FALSE), FALSE, designad, 1, FALSE)
            w <- .determine_min_variance_weight(mu = (designad@c1e + designad@c1f)/2 / sqrt(designad@n1),
                                                mu0 = 0, sigma = 1, design = designad) |> suppressMessages()
            expect_equal(.adaptively_weighted_mle(w = w, x1 = .2, x2 = .3, n1 = designad@n1,
                                                  n2 = n2_extrapol(designad, .x_to_z(x = .2, n = designad@n1, mu0 = 0, sigma = 1)), sigma = 1),
                         awsm$g2(smean1 = .2, smean2 = .3, n1 = designad@n1, n2 = n2_extrapol(designad, smean_to_z(.2, designad@n1, 1, FALSE)),
                                 sigma = 1, two_armed = FALSE),
                         tolerance = 1e-2)
          })
test_that("bias reduced estimator agrees with reference implementation.",
          {
            sw <-  get_stagewise_estimators(BiasReduced(), Normal(FALSE), FALSE, designad, 1, FALSE)
            expect_equal(.bias_reduced(mle = .3, mu = 0, mu0 = 0, sigma = 1, iterations = 1, design = designad),
                         sw$g2(designad, .3, .3, designad@n1, n2_extrapol(designad, smean_to_z(.3, designad@n1, 1, FALSE)), 1, FALSE),
                         tolerance = 1e-3)
          })
test_that("median unbiased (MLE ordering) estimator agrees with reference implementation.",
          {
            med <- get_stagewise_estimators(MedianUnbiasedMLEOrdering(), Normal(FALSE), FALSE, designad, 1, FALSE)
            expect_equal(.median_unbiased_ml(x1 = .2, x2 = .2, mu0 = 0, sigma = 1, designad),
                         med$g2(designad, .2, .2, designad@n1, n2_extrapol(designad, smean_to_z(.2, designad@n1, 1, FALSE)), 1, FALSE),
                         tolerance = 1e-3)
          })
test_that("median unbiased (LR ordering) estimator agrees with reference implementation.",
          {
            med <- get_stagewise_estimators(MedianUnbiasedLikelihoodRatioOrdering(), Normal(FALSE), FALSE, designad, 1, FALSE)
            expect_equal(.median_unbiased_lr(x1 = .2, x2 = .2, mu0 = 0, sigma = 1, designad),
                         med$g2(designad, .2, .2, designad@n1, n2_extrapol(designad, smean_to_z(.2, designad@n1, 1, FALSE)), 1, FALSE),
                         tolerance = 1e-3)
          })
test_that("median unbiased (ST ordering) estimator agrees with reference implementation.",
          {
            med <- get_stagewise_estimators(MedianUnbiasedScoreTestOrdering(), Normal(FALSE), FALSE, designad, 1, FALSE)
            expect_equal(.median_unbiased_st(x1 = .2, x2 = .2, mu0 = 0, sigma = 1, designad),
                         med$g2(designad, .2, .2, designad@n1, n2_extrapol(designad, smean_to_z(.2, designad@n1, 1, FALSE)), 1, FALSE),
                         tolerance = 1e-3)
          })
test_that("median unbiased (SWCF ordering) estimator agrees with reference implementation.",
          {
            med <- get_stagewise_estimators(MedianUnbiasedStagewiseCombinationFunctionOrdering(), Normal(FALSE), FALSE, designad, 1, FALSE)
            expect_equal(.median_unbiased_swcf(x1 = .2, x2 = .2, mu0 = 0, sigma = 1, designad),
                         med$g2(designad, .2, .2, designad@n1, n2_extrapol(designad, smean_to_z(.2, designad@n1, 1, FALSE)), 1, FALSE),
                         tolerance = 1e-3)
          })
test_that("median unbiased (Neyman Pearson ordering) estimator agrees with reference implementation.",
          {
            p <- get_stagewise_estimators(NeymanPearsonOrderingPValue(0, 0.4), Normal(FALSE), FALSE, designad, 1, FALSE)
            expect_equal(.p_np(x1 = .2, x2 = .2, mu = 0, mu0 = 0, mu1 = .4, sigma = 1, design = designad),
                         p$g2(designad, .2, .2, designad@n1, n2_extrapol(designad, smean_to_z(.2, designad@n1, 1, FALSE)),0, 1, FALSE),
                         tolerance = 1e-3)
          })
test_that("P-value (Neyman-Pearson ordering) agrees with reference implementation.",
          {
            p <- get_stagewise_estimators(NeymanPearsonOrderingPValue(0, 0.4), Normal(FALSE), FALSE, designad, 1, FALSE)
            expect_equal(.p_np(x1 = .2, x2 = .2, mu = 0, mu0 = 0, mu1 = .4, sigma = 1, design = designad),
                         p$g2(designad, .2, .2, designad@n1, n2_extrapol(designad, smean_to_z(.2, designad@n1, 1, FALSE)),0, 1, FALSE),
                         tolerance = 1e-3)
          })


test_that("CI (MLE ordering) agrees with reference implementation.",
          {
            ci <- get_stagewise_estimators(MLEOrderingCI(), Normal(FALSE), FALSE, designad, 1, FALSE)
            ref_ci <- .confidence_interval_ml(x1 = .2, x2 = .2, mu0 = 0, sigma = 1, alpha = 0.05, design=designad)
            expect_equal(ref_ci[1],
                         ci$l2(designad, .2, .2, designad@n1, n2_extrapol(designad, smean_to_z(.2, designad@n1, 1, FALSE)), 1, FALSE),
                         tolerance = 1e-2)
            expect_equal(ref_ci[2],
                         ci$u2(designad, .2, .2, designad@n1, n2_extrapol(designad, smean_to_z(.2, designad@n1, 1, FALSE)), 1, FALSE),
                         tolerance = 1e-2)
          })

test_that("CI (LR ordering) agrees with reference implementation.",
          {
            ci <- get_stagewise_estimators(LikelihoodRatioOrderingCI(), Normal(FALSE), FALSE, designad, 1, FALSE)
            ref_ci <- .confidence_interval_lr(x1 = .2, x2 = .2, mu0 = 0, sigma = 1, alpha = 0.05, design=designad)
            expect_equal(ref_ci[1],
                         ci$l2(designad, .2, .2, designad@n1, n2_extrapol(designad, smean_to_z(.2, designad@n1, 1, FALSE)), 1, FALSE),
                         tolerance = 1e-2)
            expect_equal(ref_ci[2],
                         ci$u2(designad, .2, .2, designad@n1, n2_extrapol(designad, smean_to_z(.2, designad@n1, 1, FALSE)), 1, FALSE),
                         tolerance = 1e-2)
          })
test_that("CI (ST ordering) agrees with reference implementation.",
          {
            ci <- get_stagewise_estimators(ScoreTestOrderingCI(), Normal(FALSE), FALSE, designad, 1, FALSE)
            ref_ci <- .confidence_interval_st(x1 = .2, x2 = .2, mu0 = 0, sigma = 1, alpha = 0.05, design=designad)
            expect_equal(ref_ci[1],
                         ci$l2(designad, .2, .2, designad@n1, n2_extrapol(designad, smean_to_z(.2, designad@n1, 1, FALSE)), 1, FALSE),
                         tolerance = 1e-1)
            expect_equal(ref_ci[2],
                         ci$u2(designad, .2, .2, designad@n1, n2_extrapol(designad, smean_to_z(.2, designad@n1, 1, FALSE)), 1, FALSE),
                         tolerance = 1e-1)
          })
test_that("CI (SWCF ordering) agrees with reference implementation.",
          {
            ci <- get_stagewise_estimators(StagewiseCombinationFunctionOrderingCI(), Normal(FALSE), FALSE, designad, 1, FALSE)
            ref_ci <- .confidence_interval_swcf(x1 = .2, x2 = .2, mu0 = 0, sigma = 1, alpha = 0.05, design=designad)
            expect_equal(ref_ci[1],
                         ci$l2(designad, .2, .2, designad@n1, n2_extrapol(designad, smean_to_z(.2, designad@n1, 1, FALSE)), 1, FALSE),
                         tolerance = 1e-2)
            expect_equal(ref_ci[2],
                         ci$u2(designad, .2, .2, designad@n1, n2_extrapol(designad, smean_to_z(.2, designad@n1, 1, FALSE)), 1, FALSE),
                         tolerance = 1e-2)
          })


test_that("CI (naive) agrees with reference implementation.",
          {
            ci <- get_stagewise_estimators(NaiveCI(), Normal(FALSE), FALSE, designad, 1, FALSE)
            ref_ci <- .naive_confidence_interval(x1 = .2, x2 = .2, designad@n1, n2_extrapol(designad, smean_to_z(.2, designad@n1, 1, FALSE)), sigma = 1, alpha = 0.05)
            expect_equal(ref_ci[1],
                         ci$l2(.2, .2, designad@n1, n2_extrapol(designad, smean_to_z(.2, designad@n1, 1, FALSE)), 1, FALSE),
                         tolerance = 1e-2)
            expect_equal(ref_ci[2],
                         ci$u2(.2, .2, designad@n1, n2_extrapol(designad, smean_to_z(.2, designad@n1, 1, FALSE)), 1, FALSE),
                         tolerance = 1e-2)
          })

test_that("CI (repeated) agrees with reference implementation.",
          {
            ci <- get_stagewise_estimators(RepeatedCI(), Normal(FALSE), FALSE, designad, 1, FALSE)
            ref_ci <- .repeated_confidence_interval(x1 = .2, x2 = .2, mu0 = 0, sigma = 1, design=designad)
            expect_equal(ref_ci[1],
                         ci$l2(designad, .2, .2, designad@n1, n2_extrapol(designad, smean_to_z(.2, designad@n1, 1, FALSE)), 1, FALSE),
                         tolerance = 1e-1)
            expect_equal(ref_ci[2],
                         ci$u2(designad, .2, .2, designad@n1, n2_extrapol(designad, smean_to_z(.2, designad@n1, 1, FALSE)), 1, FALSE),
                         tolerance = 1e-1)
          })

