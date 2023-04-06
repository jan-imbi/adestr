p_np <- function(design, smean, n, mu, mu0, mu1, sigma, two_armed, tol = getOption("adestr_tol_inner", default = .adestr_options[["adestr_tol_inner"]]), maxEval = getOption("adestr_maxEval_inner", default = .adestr_options[["adestr_maxEval_inner"]]), absError = getOption("adestr_absError_inner", default = .adestr_options[["adestr_absError_inner"]]), ...) {
  n1 <- n1(design, round=FALSE)
  se1 <- sigma_to_se(sigma, n1, two_armed)
  muse1 <- mu/se1
  smean1_sol <- (2L*smean*n + mu0*(n1 - n) + mu1*(n1 - n))/(2L*n1)
  z1 <- smean_to_z(smean1_sol, n1, sigma, two_armed)
  # browser()
  if (z1 < design@c1f) {
    p_part_first_stage <- pnorm(design@c1f, mean = muse1) - pnorm(z1, mean = muse1) + pnorm(design@c1e, mean = muse1, lower.tail = FALSE)
  } else if (z1 > design@c1e) {
    p_part_first_stage <- pnorm(z1, mean = muse1, lower.tail = FALSE)
  } else {
    p_part_first_stage <- pnorm(design@c1e, mean = muse1, lower.tail = FALSE)
  }
  const <- n*(2*smean - mu0 - mu1)
  p_part_second_stage <- .hcubature(f = \(x) {
    n2 <- n2_extrapol(design, x[1L,])
    se2 <- sigma * sqrt((1L + two_armed) /  n2)
    smean1_prime <- z_to_smean(x[1L,], n1, sigma, two_armed)
    nn <- (n1 + n2)
    f1_kv(x[1L,,drop=FALSE], n1, mu, sigma, two_armed) *
      pnorm( (const + mu0*nn + mu1*nn - 2*n1*smean1_prime)  / (2 * n2 * se2),
             mean = mu/se2,lower.tail = FALSE)
  },
  lowerLimit = c(design@c1f),
  upperLimit = c(design@c1e),
  tol = tol,
  maxEval = maxEval,
  absError = absError,
  vectorInterface = TRUE)$integral
  return(p_part_first_stage + p_part_second_stage)
}
p1_np <- function(design, smean1, n1, mu, mu0, mu1, sigma, two_armed, ...) p_np(design, smean1, n1, mu, mu0, mu1, sigma, two_armed, ...)
p2_np <- function(design, smean1, smean2, n1, n2, mu, mu0, mu1, sigma, two_armed, ...) p_np(design, smeans_to_smean(smean1, smean2, n1, n2), n1+n2, mu, mu0, mu1, sigma, two_armed, ...)

find_root_p1_np <- function(design, smean1, n1, mu0, mu1, sigma, two_armed, p_boundary, tol = getOption("adestr_tol_roots", default = .adestr_options[["adestr_tol_roots"]]), maxiter = getOption("adestr_maxiter_roots", default = .adestr_options[["adestr_maxiter_roots"]]), ...) {
  se1 <- sigma_to_se(sigma, n1, two_armed)
  initial_guess <- smean1 + qnorm(p_boundary, sd = se1)
  uniroot(f = \(x) p1_np(design = design, smean1 = smean1,
                         n1 = n1, mu = x, mu0 = mu0 +x, mu1 = mu1 +x, sigma = sigma, two_armed = two_armed,
                         tol = getOption("adestr_tol_inner", default = .adestr_options[["adestr_tol_inner"]]), maxEval = getOption("adestr_maxEval_inner", default = .adestr_options[["adestr_maxEval_inner"]]), absError = getOption("adestr_absError_inner", default = .adestr_options[["adestr_absError_inner"]])) - p_boundary,
          interval = c(initial_guess - se1,
                       initial_guess + se1),
          tol = tol,
          maxiter = maxiter,
          extendInt = "yes")$root
}
find_root_p2_np <- function(design, smean1, smean2, n1, n2, mu0, mu1, sigma, two_armed, p_boundary, tol = getOption("adestr_tol_roots", default = .adestr_options[["adestr_tol_roots"]]), maxiter = getOption("adestr_maxiter_roots", default = .adestr_options[["adestr_maxiter_roots"]]), ...) {
  se <- sigma_to_se(sigma, n1 + n2, two_armed)
  initial_guess <- smeans_to_smean(smean1, smean2, n1, n2) + qnorm(p_boundary, sd = se)
  uniroot(f = \(x) p2_np(design = design, smean1 = smean1, smean2 = smean2,
                         n1 = n1, n2 = n2, mu = x, mu0 = mu0 + x, mu1 = mu1 + x, sigma = sigma, two_armed = two_armed,
                         tol = getOption("adestr_tol_inner", default = .adestr_options[["adestr_tol_inner"]]), maxEval = getOption("adestr_maxEval_inner", default = .adestr_options[["adestr_maxEval_inner"]]), absError = getOption("adestr_absError_inner", default = .adestr_options[["adestr_absError_inner"]])) - p_boundary,
          interval = c(initial_guess - se,
                       initial_guess + se),
          tol = tol,
          maxiter = maxiter,
          extendInt = "yes")$root
}

setClass("NeymanPearsonOrderingPValue", contains = "VirtualPValue", slots = c(mu0 = "numeric", mu1 = "numeric"))
#' @param mu0 expected value of the normal distribution under the null hypothesis.
#' @param mu1 expected value of the normal distribution under the null hypothesis.
#' @rdname PValue-class
#' @export
NeymanPearsonOrderingPValue <- function(mu0 = 0, mu1 = 0.4) new("NeymanPearsonOrderingPValue", mu0 =  mu0, mu1 = mu1, label = paste0("Neyman-Pearson test ordering (mu0=",format(mu0),", mu1=",format(mu1),")") )
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("NeymanPearsonOrderingPValue", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            g1 <- Vectorize(\(design, smean1, n1, sigma, two_armed, ...) p1_np(design, smean1, n1, mu = 0, mu0=estimator@mu0, mu1=estimator@mu1, sigma, two_armed, ...), c("smean1"))
            g2 <- Vectorize(\(design, smean1, smean2, n1, n2, sigma, two_armed, ...) p2_np(design, smean1, smean2, n1, n2, mu = 0, mu0=estimator@mu0, mu1=estimator@mu1, sigma, two_armed, ...), c("smean1", "smean2", "n2"))
            list(g1 = g1,
                 g2 = g2)
          })
setClass("NeymanPearsonOrderingCI", contains = "VirtualIntervalEstimator", slots = c(mu0 = "numeric", mu1 = "numeric"))
#' @inheritParams NeymanPearsonOrderingPValue
#' @rdname IntervalEstimator-class
#' @export
NeymanPearsonOrderingCI <- function(two_sided = TRUE, mu0 = 0, mu1 = 0.4) new("NeymanPearsonOrderingCI", mu0 = mu0, mu1 = mu1, two_sided = two_sided, label = paste0("Neyman-Pearson test ordering (mu0=",format(mu0),", mu1=",format(mu1),")"))
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("NeymanPearsonOrderingCI", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            alpha <- adoptr_alpha_shifted_design_kv(design, 0, 0, 0)
            l1 <- Vectorize(\(design, smean1, n1, sigma, two_armed, ...) find_root_p1_np(design, smean1, n1, mu0=estimator@mu0, mu1=estimator@mu1, sigma, two_armed, p_boundary = alpha, ...), c("smean1"))
            l2 <- Vectorize(\(design, smean1, smean2, n1, n2, sigma, two_armed, ...) find_root_p2_np(design, smean1, smean2, n1, n2, mu0=estimator@mu0, mu1=estimator@mu1, sigma, two_armed, p_boundary = alpha, ...), c("smean1", "smean2", "n2"))
            if (estimator@two_sided){
              u1 <- Vectorize(\(design, smean1, n1, sigma, two_armed, ...) find_root_p1_np(design, smean1, n1, mu0=estimator@mu0, mu1=estimator@mu1, sigma, two_armed, p_boundary = 1-alpha, ...), c("smean1"))
              u2 <- Vectorize(\(design, smean1, smean2, n1, n2, sigma, two_armed, ...) find_root_p2_np(design, smean1, smean2, n1, n2, mu0=estimator@mu0, mu1=estimator@mu1, sigma, two_armed, p_boundary = 1-alpha, ...), c("smean1", "smean2", "n2"))
            } else {
              u1 <- \(design, smean1, sigma, two_armed, ...) Inf
              u2 <- \(design, smean1, smean2, sigma, two_armed, ...) Inf
            }
            return(list(
              l1 = l1,
              u1 = u1,
              l2 = l2,
              u2 = u2
            ))
          })
setClass("MedianUnbiasedNeymanPearsonOrdering", contains = "VirtualPointEstimator", slots = c(mu0 = "numeric", mu1 = "numeric"))
#' @inheritParams NeymanPearsonOrderingPValue
#' @rdname PointEstimator-class
#' @export
MedianUnbiasedNeymanPearsonOrdering <- function(mu0 = 0, mu1 = 0.4) new("MedianUnbiasedNeymanPearsonOrdering", mu0 = mu0, mu1 = mu1, label = paste0("Median unbiased (Neyman-Pearson test ordering, mu0=",format(mu0),", mu1=",format(mu1),")"))
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("MedianUnbiasedNeymanPearsonOrdering", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            list(g1 = Vectorize(\(design, smean1, n1, sigma, two_armed, ...) find_root_p1_np(design, smean1, n1, mu0=estimator@mu0, mu1=estimator@mu1, sigma, two_armed, p_boundary = 0.5, ...), c("smean1")),
                 g2 = Vectorize(\(design, smean1, smean2, n1, n2, sigma, two_armed, ...) find_root_p2_np(design, smean1, smean2, n1, n2, mu0=estimator@mu0, mu1=estimator@mu1, sigma, two_armed, p_boundary = 0.5, ...), c("smean1", "smean2", "n2")))
          })
setClass("MidpointNeymanPearsonOrderingCI", contains = "VirtualPointEstimator")
#' @inheritParams NeymanPearsonOrderingPValue
#' @rdname PointEstimator-class
#' @export
MidpointNeymanPearsonOrderingCI <- function() new("MidpointNeymanPearsonOrderingCI", label = "Midpoint of Neyman-Pearson ordering CI")
#' @rdname get_stagewise_estimators
setMethod("get_stagewise_estimators", signature("MidpointNeymanPearsonOrderingCI", "Normal"),
          function(estimator,
                   data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   sigma,
                   exact) {
            cis <- get_stagewise_estimators(estimator = NeymanPearsonOrderingCI(two_sided = TRUE),
                                            data_distribution =  data_distribution,
                                            use_full_twoarm_sampling_distribution = use_full_twoarm_sampling_distribution,
                                            design = design, sigma = sigma, exact = exact)
            l1 <- cis$l1
            l2 <- cis$l2
            u1 <- cis$u1
            u2 <- cis$u2
            g1 <- \(...) (l1(...) + u1(...))/2
            g2 <- \(...) (l2(...) + u2(...))/2
            return(list(
              g1 = g1,
              g2 = g2
            ))
          })
c2_np <- function(design, z1, mu, mu0, mu1, sigma, two_armed, alpha, tol_root = getOption("adestr_tol_roots", default = .adestr_options[["adestr_tol_roots"]]), tol = getOption("adestr_tol_inner", default = .adestr_options[["adestr_tol_inner"]]), maxEval = getOption("adestr_maxEval_inner", default = .adestr_options[["adestr_maxEval_inner"]]), absError = getOption("adestr_absError_inner", default = .adestr_options[["adestr_absError_inner"]]), ...){
  n1 <- n1(design, round=FALSE)
  se1 <- sigma_to_se(sigma, n1, two_armed)
  initialz <- qnorm(alpha/2, lower.tail = FALSE)
  z1_ord <- order(z1)
  z1 <- z1[z1_ord]
  res <- sapply(z1,
         \(z) {
           n2 <- n2_extrapol(design, z)
           n <- n1 + n2
           smean1 <- z_to_smean(z, n1, sigma, two_armed)
           root <- uniroot(\(y) {
             p_np(
               design = design,
               smean = smeans_to_smean(
                 smean1,
                 z_to_smean(y, n2, sigma, two_armed),
                 n1, n2),
               n = n, mu = mu, mu0=mu0, mu1 = mu1,
               sigma = sigma, two_armed = two_armed,
               tol = tol, maxEval = maxEval, absError = absError
               ) - alpha
           } , lower = initialz - 1, upper = initialz + 1, extendInt = "yes", tol = tol_root)$root
           initialz <<- root
           root
           })
  res <- res[order(z1_ord)]
  res
}
p_np_posterior <- function(design, smean, n, posterior_H0, posterior_H1,
                           additional_condition = \(smean_input, n_input, smean_integrand, n_integrand, ...) smean_integrand>=smean_input,
                           sigma, two_armed, tol = getOption("adestr_tol_inner", default = .adestr_options[["adestr_tol_inner"]]), maxEval = getOption("adestr_maxEval_inner", default = .adestr_options[["adestr_maxEval_inner"]]), absError = getOption("adestr_absError_inner", default = .adestr_options[["adestr_absError_inner"]]), ...) {
  logLR <- function(smean, n) log(posterior_H1(smean, n)) - log(posterior_H0(smean, n))
  n1 <- n1(design, round=FALSE)
  qInf_pH0 <- find_approx0_boundaries(posterior_H0, n=n)
  logLR_boundary <- logLR(smean, n)
  smean_cf <- z_to_smean(design@c1f, n1, sigma, two_armed)
  smean_ce <- z_to_smean(design@c1e, n1, sigma, two_armed)
  futility_p <- .hcubature(\(x) (logLR(x, n1) > logLR_boundary) *
                             additional_condition(smean_input = smean, n_input = n, smean_integrand = x, n_integrand = n1) *
                             posterior_H0(x, n1),
                           lowerLimit = c(qInf_pH0[1L]),
                           upperLimit = c(smean_cf),
                           tol = tol,
                           maxEval = maxEval,
                           absError = absError,
                           vectorInterface = TRUE)$integral
  efficacy_p <- .hcubature(\(x) (logLR(x, n1) > logLR_boundary) *
                             additional_condition(smean_input = smean, n_input = n, smean_integrand = x, n_integrand = n1) *
                             posterior_H0(x, n1),
                           lowerLimit = c(smean_ce),
                           upperLimit = c(qInf_pH0[2L]),
                           tol = tol,
                           maxEval = maxEval,
                           absError = absError,
                           vectorInterface = TRUE)$integral
  continuation_p <- .hcubature(\(x) {
    z1 <- smean_to_z(x[1L,,drop=FALSE], n1, sigma, two_armed)
    n2 <- n2_extrapol(design, z1)
    nn <- n1 + n2
    smean_integrand <- (x[1L,,drop=FALSE] * n1 + x[2L,,drop=FALSE] * n2) / nn
    (logLR(smean_integrand, nn) > logLR_boundary) *
      additional_condition(smean_input = smean, n_input = n, smean_integrand = smean_integrand, n_integrand = nn) *
      posterior_H0(x[1L,,drop=FALSE], n1) * posterior_H0(x[2L,,drop=FALSE], n2)
  },
  lowerLimit = c(smean_cf, qInf_pH0[1L]),
  upperLimit = c(smean_ce, qInf_pH0[2L]),
  tol = tol,
  maxEval = maxEval,
  absError = absError,
  vectorInterface = TRUE)$integral
  return(futility_p + efficacy_p + continuation_p)
}


# Old stuff, probably delete soon. Maybe monoLR is useful...
find_approx0_boundaries <- function(density, start = 0, tol = 1e-5, ...) {
  left <- uniroot(\(x) density(x, ...) - tol, interval = c(start-tol, start), extendInt = "upX", tol = tol)$root
  right <- uniroot(\(x) density(x, ...) - tol, interval = c(start, start+tol), extendInt = "downX", tol = tol)$root
  c(left, right)
}
find_np_min <- function(posterior_H0, posterior_H1) {
  logLR <- function(smean, n) log(posterior_H1(smean, n)) - log(posterior_H0(smean, n))
  optimize(logLR())
}
p_np_posterior_monoLR <- function(design, smean, n, posterior_H0, posterior_H1, starting_point, sigma, two_armed, tol = getOption("adestr_tol_inner", default = .adestr_options[["adestr_tol_inner"]]), maxEval = getOption("adestr_maxEval_inner", default = .adestr_options[["adestr_maxEval_inner"]]), absError = getOption("adestr_absError_inner", default = .adestr_options[["adestr_absError_inner"]]), ...) {
  logLR <- function(smean, n) log(posterior_H1(smean, n)) - log(posterior_H0(smean, n))
  n1 <- n1(design, round=FALSE)
  qInf_pH1 <- find_approx0_boundaries(posterior_H1, n=n)
  logLR_boundary <- logLR(smean, n)
  smean_cf <- z_to_smean(design@c1f, n1, sigma, two_armed)
  smean_ce <- z_to_smean(design@c1e, n1, sigma, two_armed)
  ppd_cf <- logLR(smean_cf, n1) - logLR_boundary
  ppd_ce <- logLR(smean_ce, n1) - logLR_boundary
  if (ppd_cf >= 0) {
    smean1_sol <- uniroot(
      \(x) logLR(x, n1) - logLR_boundary,
      interval = c(qInf_pH1[1L], smean_cf),
      extendInt = "upX"
    )$root
    p1 <- .hcubature(\(x) posterior_H0(x, n1),
              lowerLimit = c(smean1_sol),
              upperLimit = c(smean_cf),
              tol = tol,
              maxEval = maxEval,
              absError = absError,
              vectorInterface = TRUE)$integral
    p2 <- .hcubature(\(x) posterior_H0(x, n1),
                     lowerLimit = c(smean_ce),
                     upperLimit = c(qInf_pH1[2L]),
                     tol = tol,
                     maxEval = maxEval,
                     absError = absError,
                     vectorInterface = TRUE)$integral
    p_part_first_stage <- p1 + p2
  } else if (ppd_ce <= 0) {
    smean1_sol <- uniroot(
      \(x) logLR(x, n1) - logLR_boundary,
      interval = c(z_to_smean(design@c1e, n1, sigma, two_armed), qInf_pH1[2L]),
      extendInt = "upX"
    )$root
    p_part_first_stage <- .hcubature(\(x) posterior_H0(x, n1),
                                     lowerLimit = c(smean1_sol),
                                     upperLimit = c(qInf_pH1[2L]),
                                     tol = tol,
                                     maxEval = maxEval,
                                     absError = absError,
                                     vectorInterface = TRUE)$integral
  } else {
    # smean1_sol <- z_to_smean(design@c1e, n1, sigma, two_armed)
    p_part_first_stage <- .hcubature(\(x) posterior_H0(x, n1),
                                     lowerLimit = c(smean_ce),
                                     upperLimit = c(qInf_pH1[2L]),
                                     tol = tol,
                                     maxEval = maxEval,
                                     absError = absError,
                                     vectorInterface = TRUE)$integral
  }
  p_part_second_stage <- .hcubature(f = \(x) {
    z1 <- smean_to_z(x[1L,,drop=FALSE], n1, sigma, two_armed)
    n2 <- n2_extrapol(design, z1)
    nn <- n1 + n2
    smean2_sol <- sapply(
      seq_along(nn),
      \(i) {uniroot(
        \(yy) {logLR(smeans_to_smean(x[1L,i], yy, n1, n2[i]), nn[i]) - logLR_boundary},
        interval = c(qInf_pH1[1L], qInf_pH1[2L]),
        extendInt = "yes"
      )$root}
    )
    conditional_ps <- sapply(
      smean2_sol, \(y) {
        .hcubature(\(yy) posterior_H0(yy, n2),
          lowerLimit = y,
          upperLimit = qInf_pH1[2L],
          tol = tol,
          maxEval = maxEval,
          absError = absError,
          vectorInterface = TRUE)$integral
      })
    posterior_H0(x[1L,,drop=FALSE], n1) * conditional_ps
  },
  lowerLimit = c(smean_cf),
  upperLimit = c(smean_ce),
  tol = tol,
  maxEval = maxEval,
  absError = absError,
  vectorInterface = TRUE)$integral
  return(p_part_first_stage + p_part_second_stage)
}
check_monotone_np <- function(design, posterior_H0, posterior_H1) {
  n1 <- n1(design, round=FALSE)
  n2_candidates <- c(n1, unique(round(n2_extrapol(design, seq(design@c1f, design@c1e, length.out = 1e4)))))
  qInf_zH0 <- find_approx0_boundaries(posterior_H0, n = n1)
  # qInf_zH1 <- find_approx0_boundaries(posterior_H1, n = n1)
  qInf <- qInf_zH0
  smean_candidates <- seq(qInf[1], qInf[2], length.out=1e4)
  logLR <- function(smean, n) log(posterior_H1(smean, n)) - log(posterior_H0(smean, n))
  LRs <- sapply(
    n2_candidates,
      \(n) {
        logLR(smean_candidates, n)
      }
    )
  # LRs
  if (any(apply(LRs, 2L, \(x) x!=sort(x))))
    stop("Likelihood ratio is not monotone.")
  TRUE
}
c2_np_monoLR <- function(design, z1, posterior_H0, posterior_H1, sigma, two_armed, alpha){
  n1 <- n1(design, round=FALSE)
  qInf_zH0 <- smean_to_z(find_approx0_boundaries(posterior_H0), n1, sigma, two_armed)
  sapply(z1,
         \(z) {
           n2 <- n2_extrapol(design, z)
           n <- n1 + n2
           smean1 <- z_to_smean(z, n1, sigma, two_armed)
           uniroot(\(y) {
           p_np_posterior_monoLR(
             design = design,
             smean = smeans_to_smean(
               smean1,
               z_to_smean(y, n2, sigma, two_armed),
               n1, n2),
             n = n, posterior_H0 = posterior_H0, posterior_H1 = posterior_H1,
             sigma = sigma, two_armed = two_armed) - alpha
         } , lower = qInf_zH0[1], upper = qInf_zH0[2], extendInt = "yes", tol = 1e-6)$root}
  )
}

p_np_shifted <- function(design, smean, n, mu, mu0, mu1, sigma, two_armed, tol = getOption("adestr_tol_inner", default = .adestr_options[["adestr_tol_inner"]]), maxEval = getOption("adestr_maxEval_inner", default = .adestr_options[["adestr_maxEval_inner"]]), absError = getOption("adestr_absError_inner", default = .adestr_options[["adestr_absError_inner"]]), ...) {
  n1 <- n1(design, round=FALSE)
  se1 <- sigma_to_se(sigma, n1, two_armed)
  z1 <- smean_to_z(smean, n1, sigma, two_armed)
  muse1 <- mu/se1
  smean1_sol <- (2L*smean*n + mu0*(n1 - n) + mu1*(n1 - n))/(2L*n1)
  z1_sol <- smean_to_z(smean1_sol, n1, sigma, two_armed)

  if (z1 < design@c1f) {
    p_part_first_stage <- pnorm(z1, mean=muse1, lower.tail = FALSE)
    p_part_second_stage <- 0
  } else {
    # if (z1 > design@c1e) {
      p_part_first_stage <- pnorm(z1, mean = muse1, lower.tail = FALSE)
    # } else {
    #   p_part_first_stage <- pnorm(design@c1e, mean = muse1, lower.tail = FALSE)
    # }
    const <- n*(2*smean - mu0 - mu1)
    # const <- n * (2 * smean  - mu1)
    p_part_second_stage <- .hcubature(f = \(x) {
      n2 <- n2_extrapol(design, x[1L,])
      se2 <- sigma * sqrt((1L + two_armed) /  n2)
      smean1_prime <- z_to_smean(x[1L,], n1, sigma, two_armed)
      nn <- (n1 + n2)
      f1_kv(x[1L,,drop=FALSE], n1, mu, sigma, two_armed) *
        pnorm( (const + mu0*nn + mu1*nn - 2*n1*smean1_prime)  / (2 * n2 * se2),
               # pnorm( (const + mu1 * (n1 + n2)  - 2 * n1 * smean1_prime)  / (2 * n2 * se2),
               mean = mu/se2,lower.tail = FALSE)
    },
    lowerLimit = c(design@c1f),
    upperLimit = c(min(design@c1e, z1)),
    tol = tol,
    maxEval = maxEval,
    absError = absError,
    vectorInterface = TRUE)$integral
  }
  return(p_part_first_stage + p_part_second_stage)
}

