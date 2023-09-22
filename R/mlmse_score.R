### Uncomment this once adoptr is back on CRAN ###

# setClass("MLMSE",
#          slots = c(
#            mu = "numeric",
#            sigma = "numeric",
#            two_armed = "logical",
#            tol = "numeric",
#            maxEval = "numeric",
#            absError = "numeric"
#          ),
#          contains = "UnconditionalScore")
# MLMSE <- function(mu = 0, sigma = 1, two_armed = FALSE, tol = 1e-5, maxEval = 1e7, absError = 1e-7) {
#   new("MLMSE",
#       mu = mu,
#       sigma = sigma,
#       two_armed = two_armed,
#       tol = tol,
#       maxEval = maxEval,
#       absError = absError,
#       label = "E[|mu_ML - mu|^2]")
# }
# setMethod("evaluate",
#           signature("MLMSE", "TwoStageDesign"),
#           function(s, design, ...) {
#             mean(
#               evaluate_estimator(
#                 score = MSE(),
#                 estimator = SampleMean(),
#                 data_distribution = Normal(two_armed = s@two_armed),
#                 design = design,
#                 mu = s@mu,
#                 sigma = s@sigma,
#                 tol = s@tol,
#                 maxEval = s@maxEval
#               )@results$MSE
#             )
#           })
#
#
