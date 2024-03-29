get_norm_inf <- function(infq, means, left=TRUE) {
  extremum <- if (left) min else max
  extremum(qnorm(infq, mean = means, lower.tail = left))
}
get_t_inf <- function(infq, dfs, means, left=TRUE) {
  extremum <- if (left) min else max
  extremum(suppressWarnings(qt(infq, dfs, means, lower.tail = left)))
}
get_chi_inf <- function(infq, dfs, variance, left=TRUE) {
  extremum <- if (left) min else max
  extremum(qchisq(infq, dfs, lower.tail = left) * variance / dfs)
}

# Other stuff
z_to_smean <- function(z, n, sigma, two_armed) z * sigma * sqrt((1L + two_armed) / n)
t_to_smean <- function(t, svar, n, two_armed) t * sqrt(svar * (1L + two_armed) / n)
t_to_z <- function(t, svar, sigma, n) t * sqrt(svar)/sigma
z_to_t <- function(z, svar, sigma, n) z * sigma / sqrt(svar)
smean_to_z <- function(smean, n, sigma, two_armed) smean * sqrt(n / (1L + two_armed)) / sigma
smean_to_t <- function(smean, svar, n, two_armed) smean * sqrt(n / (svar * (1L + two_armed)))
sigma_to_se <- function(sigma, n, two_armed) sigma * sqrt((1L + two_armed) / n)
svar_to_se <- function(svar, n, two_armed) sqrt(svar * (1L + two_armed) / n)
smeans_to_smean <- function(smean1, smean2, n1, n2) (n1 * smean1 + n2 * smean2) / (n1 + n2)

z_test <- function(smean, n, sigma, two_sample) {
 if (two_sample && length(n)!=2L)
   stop("Two-sample Z-test requires that n has two entries.")
  if (two_sample)
    return(smean / (sigma * sqrt(1/n[1L] + 1/n[2L])))
  else
    return(smean / (sigma * sqrt(1/n[1L])))
}
t_test <- function(smean, svar, n, two_sample) {
  if (two_sample && (length(n)!=2L))
    stop("Two-sample t-test requires that n has two entries.")
  if (two_sample)
    return(smean / sqrt(svar * (1/n[1L] + 1/n[2L])) )
  else
    return(smean / sqrt(svar * 1/n[1L]))
}

pool_svar1_svar2 <- function(svar1, svar2, n1, n2, two_armed) {
  df1 <- (1L + two_armed) * (n1 - 1L)
  df2 <- (1L + two_armed) * (n2 - 1L)
  (svar1 * df1 + svar2 * df2)/(df1 + df2)
}
get_overall_svar_onearm <- function(smean1, svar1, smean2, svar2, n1, n2){
  n <- n1 + n2

  ssum1 <- smean1 * n1
  ssum2 <- smean2 * n2
  ssum <- ssum1 + ssum2

  smean <- (ssum1 + ssum2) / (n1 + n2)

  smeansq  <- smean^2L
  smeansq1 <- smean1^2L
  smeansq2 <- smean2^2L

  ssqsum1 <- (n1-1L) * svar1 + n1 *  smeansq1
  ssqsum2 <- (n2-1L) * svar2 + n2 *  smeansq2

  ssqsum  <- ssqsum1 + ssqsum2

  svar <- (ssqsum  - smeansq * n)  / (n1 + n2-1L)

  return(svar)
}
get_overall_svar_twoarm <- function(smean1, smean1T, svar1, smean2, smean2T, svar2, n1, n2){
  smean1C <- -smean1 +smean1T
  smean2C <- -smean2 +smean2T
  n <- n1 + n2
  ssum1T <- smean1T * n1
  ssum1C <- smean1C * n1
  ssum2T <- smean2T * n2
  ssum2C <- smean2C * n2
  ssumT <- ssum1T + ssum2T
  ssumC <- ssum1C + ssum2C
  smeanT <- ssumT / (n1 + n2)
  smeanC <- ssumC / (n1 + n2)
  smeansqT  <- smeanT^2L
  smeansq1T <- smean1T^2L
  smeansq2T <- smean2T^2L
  smeansqC  <- smeanC^2L
  smeansq1C <- smean1C^2L
  smeansq2C <- smean2C^2L
  smeansq <- smeansqT + smeansqC

  ssqsum1 <- 2*(n1-1L) * svar1 + n1 *  (smeansq1T + smeansq1C)
  ssqsum2 <- 2*(n2-1L) * svar2 + n2 *  (smeansq2T + smeansq2C)

  ssqsum  <- ssqsum1 + ssqsum2
  svar <- (ssqsum  - smeansq * n)  / (2*(n1 + n2-1L))
  return(svar)
}

#' Generate an exemplary adaptive design
#'
#' The design was optimized to minimize the expected sample size
#' under the alternative hypothesis for a one-armed trial.
#' The boundaries are chosen to control the type I error at 0.025
#' for a normally distributed test statistic (i.e. known variance).
#' For an alternative hypothesis of mu=0.4, the overall power is 80\%.
#'
#' @param two_armed (logical) determins whether the design is for one- or
#' two-armed trials.
#' @param label (optional) label to be assigned to the design.
#'
#' @returns an exemplary design of class \code{TwoStageDesign}. This object
#' contains information about the sample size recalculation rule \code{n2}, the
#' futility and efficacy boundaries \code{c1f} and \code{c1e} and the
#' second-stage rejection boundary \code{c2}.
#'
#' @export
#'
#' @examples
#' get_example_design()
#'
get_example_design <- function(two_armed = FALSE, label = NULL) {
  if (two_armed) {
    d <- TwoStageDesign(
      n1 = 56.33739084822602904978,
      c1f = 0.7907356135206976555097,
      c1e = 2.291313615804561720779,
      n2_pivots = c(
        78.74984462914770233510,
        74.46976145811771630179,
        66.54436357139142899086,
        55.54462857388867291775,
        42.83473938241603207189,
        30.34059894227031151104,
        20.48959489543554823854
      ),
      c2_pivots = c(
        2.16905734410577055726321,
        2.02492349183474207308109,
        1.77298374394898416994693,
        1.42521477360223225439029,
        0.99909735679793421070372,
        0.52314051699418129270924,
        0.07058637352917693230658
      ),
      7
    )
  } else {
    d <- TwoStageDesign(
      n1 = 28.16834031633078083701,
      c1f = 0.7907304707554818623549,
      c1e = 2.291260947864900643367,
      n2_pivots = c(
        39.39294353955478555918,
        37.23397813905835818105,
        33.27173714438612961430,
        27.77227568901122012335,
        21.41776450755991234587,
        15.17163280081247300757,
        10.25508398663193787570
      ),
      c2_pivots = c(
        2.16914648055318837194250,
        2.02493357331804890719695,
        1.77299079624771049878973,
        1.42524439642541422834654,
        0.99916431580845133098023,
        0.52325801518650127963639,
        0.07133753446126563091401
      ),
      7
    )
  }
  d <- TwoStageDesignWithCache(d)
  attr(d, "label") <- label
  d
}
#' Generate a list of estimators and p-values to use in examples
#'
#' This function generates a list of objects of class \code{\link{PointEstimator}},
#' \code{\link{IntervalEstimator}}s, and \code{\link{PValue}}s to use in
#' examples of the \code{\link{analyze}} function.
#'
#' @details
#' ## Point estimators
#' The following \code{\link{PointEstimator}}s are included:
#'  * \code{\link{SampleMean}}
#'  * \code{\link{PseudoRaoBlackwell}}
#'  * \code{\link{MedianUnbiasedLikelihoodRatioOrdering}}
#'  * \code{\link{BiasReduced}}
#'
#' ## Confidence intervals
#' The following \code{\link{IntervalEstimator}}s are included:
#'  * \code{\link{StagewiseCombinationFunctionOrderingCI}}
#'  * \code{\link{LikelihoodRatioOrderingCI}}
#'
#' ## P-Values
#' The following \code{\link{PValue}}s are included:
#'  * \code{\link{StagewiseCombinationFunctionOrderingPValue}}
#'  * \code{\link{LikelihoodRatioOrderingPValue}}
#' @md
#'
#' @param point_estimators logical indicating whether point estimators should be included in output list
#' @param interval_estimators logical indicating whether interval estimators should be included in output list
#' @param p_values logical indicating whether p-values should be included in output list
#'
#' @returns a list of \code{\link{PointEstimator}}s, \code{\link{IntervalEstimator}}s and
#' \code{\link{PValue}}.
#' @export
#'
#' @inherit analyze examples
get_example_statistics <- function(point_estimators = TRUE,
                                   interval_estimators = TRUE,
                                   p_values = TRUE) {
  point <- list(SampleMean(), PseudoRaoBlackwell(), MedianUnbiasedLikelihoodRatioOrdering(), BiasReduced())
  interval <- list(StagewiseCombinationFunctionOrderingCI(), LikelihoodRatioOrderingCI())
  p <- list(StagewiseCombinationFunctionOrderingPValue(), LikelihoodRatioOrderingPValue())
  ret <- list()
  if (point_estimators)
    ret <- c(ret, point)
  if (interval_estimators)
    ret <- c(ret, interval)
  if (p_values)
    ret <- c(ret, p)
  return(ret)
}

#' Generate the list of estimators and p-values that were used in the paper
#'
#' @param point_estimators logical indicating whether point estimators should be included in output list
#' @param interval_estimators logical indicating whether interval estimators should be included in output list
#' @param p_values logical indicating whether p-values should be included in output list
#'
#' @returns a list of \code{\link{PointEstimator}}s, \code{\link{IntervalEstimator}}s and
#' \code{\link{PValue}}.
#'
#' @inherit analyze examples
get_statistics_from_paper <- function(point_estimators = TRUE,
                                   interval_estimators = TRUE,
                                   p_values = TRUE) {
  point <- list(
    SampleMean(),
    FirstStageSampleMean(),
    MinimizePeakVariance(),
    PseudoRaoBlackwell(),
    BiasReduced(),
    MedianUnbiasedMLEOrdering(),
    MedianUnbiasedLikelihoodRatioOrdering(),
    MedianUnbiasedScoreTestOrdering(),
    MedianUnbiasedStagewiseCombinationFunctionOrdering()
  )
  interval <- list(
    NaiveCI(),
    MLEOrderingCI(),
    LikelihoodRatioOrderingCI(),
    ScoreTestOrderingCI(),
    StagewiseCombinationFunctionOrderingCI(),
    RepeatedCI()
  )
  p <- list(
    MLEOrderingPValue(),
    ScoreTestOrderingPValue(),
    #NeymanPearsonOrderingPValue(),
    StagewiseCombinationFunctionOrderingPValue(),
    LikelihoodRatioOrderingPValue())
  ret <- list()
  if (point_estimators)
    ret <- c(ret, point)
  if (interval_estimators)
    ret <- c(ret, interval)
  if (p_values)
    ret <- c(ret, p)
  return(ret)
}
