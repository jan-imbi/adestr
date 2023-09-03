design <- get_example_design()
n2_extrapol(design, 1.4)
n2(design, 1.4, round = FALSE)
integrate(.f1, lower = -Inf, upper = Inf, n1 = 30, mu = 0, mu0 = 0, sigma = 1)
hcubature(
  \(x, n1, n2, mu, mu0, sigma, design) {
    .f2(x[1,,drop=FALSE], x[2,,drop=FALSE],
        n1 = n1, n2 = n2, mu = mu, mu0 = mu0, sigma = sigma)
    },
  vectorInterface = TRUE,
  lowerLimit = c(-Inf, -Inf),
  upperLimit = c(Inf, Inf),
  n1 = 30,
  n2 = 30,
  mu = 0,
  mu0 = 0,
  sigma = 1,
  design = design
)$integral


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









