library(adestr)
designad <- get_example_design()
designgs <- readRDS("data/designgs.rds")

options(
  list(
    # Root finding inside estimators
    adestr_tol_roots = 1e-3,
    adestr_maxiter_roots = 1e3,
    # Integrals used inside estimators
    adestr_tol_inner = 1e-2,
    adestr_maxEval_inner = 1e3,
    adestr_absError_inner = 1e-5,
    # Integrals to evaluate estimators
    adestr_tol_outer = 1e-3,
    adestr_maxEval_outer = 1e4,
    adestr_absError_outer = 1e-8
  )
)
a <- evaluate_estimator(
  score = TestAgreement(),
  estimator = NaiveCI(),
  data_distribution = Normal(FALSE),
  use_full_twoarm_sampling_distribution = FALSE,
  design = designgs,
  mu = 0.15,
  sigma = 1
)
a
a@integrals


a <- evaluate_estimator(
  score = TestAgreement(),
  estimator = MLEOrderingCI(),
  data_distribution = Student(two_armed = TRUE),
  design = designad,
  mu = .3,
  sigma = 1
)
a



