# Set tolerances for unit testing
options(
  list(
    # Root finding inside estimators
    adestr_tol_roots = 1e-3,
    adestr_maxiter_roots = 1e3,
    # Integrals used inside estimators
    adestr_tol_inner = 5e-3,
    adestr_maxEval_inner = 1e3,
    adestr_absError_inner = 1e-5,
    # Integrals to evaluate estimators
    adestr_tol_outer = 5e-6,
    adestr_maxEval_outer = 3e4,
    adestr_absError_outer = 1e-8
  )
)

evaluate_estimator(MSE(), SampleMean(), Normal(FALSE), design = designad, mu=0.3, sigma = 1)
evaluate_estimator(MSE(), MedianUnbiasedMLEOrdering(), Normal(FALSE), design = designad, mu=0.3, sigma = 1)
evaluate_estimator(MSE(), BiasReduced(), Normal(FALSE), design = designad, mu=0.3, sigma = 1)

evaluate_estimator(MSE(), SampleMean(), Normal(TRUE), design = designad, mu=0.3, sigma = 1)
evaluate_estimator(MSE(), MedianUnbiasedMLEOrdering(), Normal(TRUE), design = designad, mu=0.3, sigma = 1)
evaluate_estimator(MSE(), BiasReduced(), Normal(TRUE), design = designad, mu=0.3, sigma = 1)


evaluate_estimator(MSE(), SampleMean(), Student(FALSE), design = designad, mu=0.3, sigma = 1)
evaluate_estimator(MSE(), MedianUnbiasedMLEOrdering(), Student(FALSE), design = designad, mu=0.3, sigma = 1)
evaluate_estimator(MSE(), BiasReduced(), Student(FALSE), design = designad, mu=0.3, sigma = 1)
