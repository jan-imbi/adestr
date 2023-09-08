# Author: Jan Meis
# Source: https://github.com/jan-imbi/adestr/

# This script performs the necessary calculations to reproduce
# the results from the paper (by default with lower accuracy
# save computation time).
# If you do not want to run the calculations yourself, you can get the data from
# https://github.com/jan-imbi/adestr/tree/master/data

# First, we will need to install adestr and a special version of adoptr
library(remotes)
install_github("imbi-heidelberg/adoptr", ref = "adestr_paper")
install_github("jan-imbi/adestr")
library(adestr)

# Next, we will enable parallel execution
library(future)

# Enable progress bars.
library(progressr)
handlers(global = TRUE)
handlers("progress")

## If you want to reproduce calculations on your local machine
## with lower accuracy, do this
plan(multisession, workers = parallelly::availableCores()-1L)
options(list(
  # Root finding inside estimators
  adestr_tol_roots = 1e-3,
  adestr_maxiter_roots = 1e3,
  # Integrals used inside estimators
  adestr_tol_inner = 5e-2,
  adestr_maxEval_inner = 5e2,
  adestr_absError_inner = 1e-5,
  # Integrals to evaluate estimators
  adestr_tol_outer = 1e-4,
  adestr_maxEval_outer = 3e3,
  adestr_absError_outer = 1e-6
))
## Otherwise, if you have a server and want high-accuracy result,
## do this:
# plan(multisession, workers = 45)
# options(list(
#   # Root finding inside estimators
#   adestr_tol_roots = 1e-3,
#   adestr_maxiter_roots = 1e3,
#   # Integrals used inside estimators
#   adestr_tol_inner = 5e-3,
#   adestr_maxEval_inner = 1e4,
#   adestr_absError_inner = 1e-5,
#   # Integrals to evaluate estimators
#   adestr_tol_outer = 5e-6,
#   adestr_maxEval_outer = 3e5,
#   adestr_absError_outer = 1e-8
# ))

# Calculate the optimal design parameters with adoptr
## Main paper: One-arm design assuming known variance
dir.create("data")
datadist <- Normal(two_armed = FALSE)
H_0 <- PointMassPrior(.0, 1)
H_1 <- PointMassPrior(.4, 1)
H_01 <- PointMassPrior(c(0, .4), c(.5, .5))
ess_H0 <- ExpectedSampleSize(datadist, H_0)
ess_H1 <- ExpectedSampleSize(datadist, H_1)
toer  <- Power(datadist, H_0)
power <- Power(datadist, H_1)
cp <- ConditionalPower(datadist, H_1)
initial_gs <- get_initial_design(
  theta = .4,
  alpha = .025,
  beta  = .2,
  type_design  = "group-sequential",
  dist  = datadist
)
initial_ad <- get_initial_design(
  theta = .4,
  alpha = .025,
  beta  = .2,
  type_design  = "two-stage",
  dist  = datadist
)
optad <- minimize(
  ess_H1,
  subject_to(
    power >= 0.8,
    toer  <= .025
  ),
  initial_ad,
  opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-05, maxeval = 40000)
)
designad <- adestr:::cache_design_splines(optad$design)
attr(designad, "label") <- "Adaptive"

optgs <- minimize(
  ess_H1,
  subject_to(
    power >= 0.8,
    toer  <= .025
  ),
  initial_gs,
  opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-05, maxeval = 20000)
)
designgs <- adestr:::cache_design_splines(optgs$design)
attr(designgs, "label") <- "Group-sequential"
saveRDS(designad, "data/designad.rds")
saveRDS(designgs, "data/designgs.rds")

## Supplement: Two-arm design assuming unknown variance (T-distribution)
datadist_t <- Student(two_armed = TRUE)
ess_H0_t <- ExpectedSampleSize(datadist_t, H_0)
ess_H1_t <- ExpectedSampleSize(datadist_t, H_1)
toer_t  <- Power(datadist_t, H_0)
power_t <- Power(datadist_t, H_1)
initial_gs_t <- get_initial_design(
  theta = .4,
  alpha = .025,
  beta  = .2,
  type_design  = "group-sequential",
  dist  = datadist_t
)
initial_ad_t <- get_initial_design(
  theta = .4,
  alpha = .025,
  beta  = .2,
  type_design  = "two-stage",
  dist  = datadist_t
)
optad_t <- minimize(
  ess_H1_t,
  subject_to(
    power_t >= 0.8,
    toer_t  <= .025
  ),
  initial_ad_t,
  opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-06, maxeval = 40000)
)
designad_t <- adestr:::cache_design_splines(optad_t$design)
attr(designad_t, "label") <- "Adaptive"

optgs_t <- minimize(
  ess_H1_t,
  subject_to(
    power_t >= 0.8,
    toer_t  <= .025
  ),
  initial_gs_t,
  opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-06, maxeval = 20000)
)
designgs_t <- adestr:::cache_design_splines(optgs_t$design)
attr(designgs_t, "label") <- "Group-sequential"
saveRDS(designad_t, "data/designad_t.rds")
saveRDS(designgs_t, "data/designgs_t.rds")

# Define the parameter lists for the scenarios to be evaluated
# (which estimators, assumptions on mu, which performance measures, etc.)
muvec <- sort(c(seq(-0.75, 1.32, .03),
                mean(c(designad@c1f, designad@c1e))/sqrt(designad@n1),
                mean(c(designgs@c1f, designgs@c1e))/sqrt(designgs@n1)))
saveRDS(muvec, "data/muvec.rds")
point_estimator_score_list <- list(MSE(), OverestimationProbability())
interval_estimator_score_list <- list(Coverage(), Width(), TestAgreement())
point_estimator_list <- list(SampleMean(), FirstStageSampleMean(), PseudoRaoBlackwell(), BiasReduced(), MinimizePeakVariance(),
                             MedianUnbiasedStagewiseCombinationFunctionOrdering(), MedianUnbiasedMLEOrdering(), MedianUnbiasedLikelihoodRatioOrdering(),
                             MedianUnbiasedScoreTestOrdering())
interval_estimator_list <- list(StagewiseCombinationFunctionOrderingCI(), MLEOrderingCI(), LikelihoodRatioOrderingCI(),
                                ScoreTestOrderingCI(), NeymanPearsonOrderingCI(), RepeatedCI(), NaiveCI())
dist <- list(Normal(two_armed = FALSE))
dist_t <- list(Student(two_armed = TRUE))
designs <- list(designad, designgs)
designs_t <- list(designad_t, designgs_t)

# Main paper
tab_list_main <- evaluate_scenarios_parallel(
  score_lists = list(point_estimator_score_list, interval_estimator_score_list),
  estimator_lists = list(point_estimator_list, interval_estimator_list),
  data_distribution_lists = list(dist, dist),
  design_lists = list(designs, designs),
  mu_lists = list(muvec, muvec),
  sigma_lists = list(1, 1)
)
saveRDS(tab_list_main[[1L]][,-length(tab_list_main[[1L]])],  "data/tab_est.rds")
## Somehow, one integral failed to converge properly, this fixed that
fix_bad_convergence <-
  evaluate_estimator(
    score = TestAgreement(),
    estimator = NaiveCI(),
    data_distribution = Normal(FALSE),
    use_full_twoarm_sampling_distribution = FALSE,
    design = designgs,
    mu = 0.15,
    sigma = 1,
    tol = 1e-8
  )
tmp <- tab_list_main[[2L]][,-length(tab_list_main[[2L]])]
tmp[abs(tmp$mu - 0.15) < .0001 &
      tmp$design == "Group-sequential" & tmp$estimator == "Naive CI",
    c(
      "Agreement with test decision",
      "error_Agreement with test decision",
      "functionEvaluations_Agreement with test decision"
    )] <- c(fix_bad_convergence@results$`Agreement with test decision`,
            fix_bad_convergence@integrals$`Agreement with test decision`$overall_integral$error,
            fix_bad_convergence@integrals$`Agreement with test decision`$overall_integral$functionEvaluations)
saveRDS(tmp, "data/tab_ci.rds")
rm(tmp)

# Supplement
tab_list_supp <- evaluate_scenarios_parallel(
  score_lists = list(point_estimator_score_list, interval_estimator_score_list),
  estimator_lists = list(point_estimator_list, interval_estimator_list),
  data_distribution_lists = list(dist_t, dist_t),
  design_lists = list(designs_t, designs_t),
  mu_lists = list(muvec, muvec),
  sigma_lists = list(1, 1)
)
saveRDS(tab_list_supp[[1L]][,-length(tab_list_supp[[1L]])],  "data/tab_est_t.rds")
saveRDS(tab_list_supp[[2L]][,-length(tab_list_supp[[2L]])],  "data/tab_ci_t.rds")

