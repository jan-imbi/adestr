.adestr_options <- list(
  # Root finding inside estimators
  adestr_tol_roots = eval(formals(uniroot)$tol),
  adestr_maxiter_roots = eval(formals(uniroot)$maxiter),
  # Integrals used inside estimators
  adestr_tol_inner = eval(formals(hcubature)$tol),
  adestr_maxEval_inner = eval(formals(hcubature)$maxEval),
  adestr_absError_inner = eval(formals(hcubature)$absError),
  # Integrals to evaluate estimators
  adestr_tol_outer = eval(formals(hcubature)$tol),
  adestr_maxEval_outer = eval(formals(hcubature)$maxEval),
  adestr_absError_outer = eval(formals(hcubature)$absError)
)
