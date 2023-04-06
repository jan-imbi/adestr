#' @importFrom cubature hcubature
.hcubature <- function(f,
                       lowerLimit,
                       upperLimit,
                       ...,
                       tol = 1e-05,
                       fDim = 1,
                       maxEval = 0,
                       absError = 0,
                       doChecking = FALSE,
                       vectorInterface = FALSE,
                       norm = c("INDIVIDUAL", "PAIRED", "L2", "L1", "LINF")) {
  if (any(upperLimit <= lowerLimit)) {
    return(list(
      integral = 0,
      error = 0,
      functionEvaluations = 0,
      returnCode = 0
    ))
  }
  else {
    return(
      hcubature(f = f,
              lowerLimit = lowerLimit,
              upperLimit = upperLimit,
              ...,
              tol = tol,
              fDim = fDim,
              maxEval = maxEval,
              absError = absError,
              doChecking = doChecking,
              vectorInterface = vectorInterface,
              norm = norm)
    )
  }
}
