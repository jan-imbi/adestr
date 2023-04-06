library(microbenchmark)
designad <- readRDS("F:/Documents/code_projects/adestr/data/designad.rds")

microbenchmark(
  evaluate_estimator(SoftCoverage(1e4), LikelihoodRatioOrderingCI(), Normal(FALSE), FALSE, designad, mu=0, sigma = 1),
  evaluate_estimator(Coverage(), LikelihoodRatioOrderingCI(), Normal(FALSE), FALSE, designad, mu=0, sigma = 1),
 times = 1
)


