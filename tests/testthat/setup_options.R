# Set tolerances for unit testing
options(
  list(
    # Root finding inside estimators
    adestr_tol_roots = 1e-3,
    adestr_maxiter_roots = 1e3,
    # Integrals used inside estimators
    adestr_tol_inner = 5e-3,
    adestr_maxEval_inner = 1e4,
    adestr_absError_inner = 1e-5,
    # Integrals to evaluate estimators
    adestr_tol_outer = 5e-6,
    adestr_maxEval_outer = 3e5,
    adestr_absError_outer = 1e-8
  )
)
