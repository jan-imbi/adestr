test_that("densities agree with reference implementation.",
{
  expect_equal(.f1(z1 = 1.5, n1 = 30, mu = 0, mu0 = 0, sigma = 1),
               f1_kv(z1 = 1.5, n1 = 30, mu = 0, sigma = 1, two_armed = FALSE))
  expect_equal(.f2(z1 = 1.5, z2 = 1, n1 = 30, n2 = 45, mu = 0, mu0 = 0, sigma = 1),
               f2_kv(z1 = 1.5, z2 = 1, n1 = 30, n2 = 45, mu = 0, sigma = 1, two_armed = FALSE))

})





.mle_pdf(0.5, 0, 0, 1, design)

dsmean(
  Normal(two_armed = FALSE),
  design,
  .5,
  0,
  1,
  exact = FALSE,
  combine_components = FALSE
)

.ml_expectation_derivative(0, 0, 1, design)



sw <-  get_stagewise_estimators(BiasReduced(), Normal(FALSE), FALSE, design, 1, FALSE)
sw$g2(design, .0, .0, 30, 30, 1, FALSE)
.bias_reduced(.0, 0, 0, 1, 1, design)

pml <- get_stagewise_estimators(MLEOrderingPValue(), Normal(FALSE), FALSE, design, 1, FALSE)

pml$g2(design, .2, .5, design@n1, n2(design, smean_to_z(.2, design@n1, 1, FALSE), round = FALSE), 1, FALSE)
.p_ml(.2, .5, 0, 0, 1, design)


prb <- get_stagewise_estimators(PseudoRaoBlackwell(), Normal(FALSE), FALSE, design, 1, FALSE)

prb$g2(design, .4, .4, design@n1, n2_extrapol(design, smean_to_z(.4, design@n1, 1, FALSE)),
       0, 1, FALSE)
.pseudo_rao_blackwell(.4, .4, .0, 1, design)



awsm <- get_stagewise_estimators(MinimizePeakVariance(), Normal(FALSE), FALSE, design, 1, FALSE)






