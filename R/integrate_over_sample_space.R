setGeneric("integrate_over_sample_space",
           function(data_distribution,
                    use_full_twoarm_sampling_distribution = FALSE,
                    design,
                    g1 = \(...) 1L,
                    g2 = \(...) 1L,
                    mu,
                    sigma,
                    tol = getOption("adestr_tol_outer", default = .adestr_options[["adestr_tol_outer"]]),
                    maxEval = getOption("adestr_maxEval_outer", default = .adestr_options[["adestr_maxEval_outer"]]),
                    absError = getOption("adestr_absError_outer", default = .adestr_options[["adestr_absError_outer"]]),
                    exact = FALSE,
                    early_futility_part = TRUE,
                    continuation_part = TRUE,
                    early_efficacy_part = TRUE,
                    conditional_integral = FALSE) standardGeneric("integrate_over_sample_space"))

setMethod("integrate_over_sample_space", signature = signature("Normal"),
          function(data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   g1,
                   g2,
                   mu,
                   sigma,
                   tol = getOption("adestr_tol_outer", default = .adestr_options[["adestr_tol_outer"]]),
                   maxEval = getOption("adestr_maxEval_outer", default = .adestr_options[["adestr_maxEval_outer"]]),
                   absError = getOption("adestr_absError_outer", default = .adestr_options[["adestr_absError_outer"]]),
                   exact,
                   early_futility_part,
                   continuation_part,
                   early_efficacy_part,
                   conditional_integral) {
            if (use_full_twoarm_sampling_distribution & !(data_distribution@two_armed))
              stop("data_distribution is one-armed but use_full_twoarm_sampling_distribution is true.")
            if (!use_full_twoarm_sampling_distribution) {
              return(
                int_kv(design,
                             g1,
                             g2,
                             mu,
                             sigma,
                             data_distribution@two_armed,
                             tol,
                             maxEval,
                             absError,
                             exact,
                             early_futility_part,
                             continuation_part,
                             early_efficacy_part,
                             conditional_integral))
            } else {
              return(
                int_kv_full(design,
                             g1,
                             g2,
                             mu,
                             sigma,
                             tol,
                             maxEval,
                             absError,
                             exact,
                             early_futility_part,
                             continuation_part,
                             early_efficacy_part,
                             conditional_integral))
            }
          })
setMethod("integrate_over_sample_space", signature = signature("Student"),
          function(data_distribution,
                   use_full_twoarm_sampling_distribution,
                   design,
                   g1,
                   g2,
                   mu,
                   sigma,
                   tol,
                   maxEval,
                   absError,
                   exact,
                   early_futility_part,
                   continuation_part,
                   early_efficacy_part,
                   conditional_integral) {
            if (use_full_twoarm_sampling_distribution & !(data_distribution@two_armed))
              stop("data_distribution is one-armed but use_full_twoarm_sampling_distribution is true.")
            if (!use_full_twoarm_sampling_distribution) {
              int_uv(design,
                           g1,
                           g2,
                           mu,
                           sigma,
                           data_distribution@two_armed,
                           tol,
                           maxEval,
                           absError,
                           exact,
                           early_futility_part,
                           continuation_part,
                           early_efficacy_part,
                           conditional_integral)
            } else {
              int_uv_full(design,
                          g1,
                          g2,
                          mu,
                          sigma,
                          tol,
                          maxEval,
                          absError,
                          exact,
                          early_futility_part,
                          continuation_part,
                          early_efficacy_part,
                          conditional_integral)
            }
          })

# integrate distribution of statistic for 1 unknown parameter (mu) over adoptr design
#' @importFrom cubature hcubature
int_kv <- function(design,
                   g1 = \(...) 1L,
                   g2 = \(...) 1L,
                   mu,
                   sigma,
                   two_armed,
                   tol = getOption("adestr_tol_outer", default = .adestr_options[["adestr_tol_outer"]]),
                   maxEval = getOption("adestr_maxEval_outer", default = .adestr_options[["adestr_maxEval_outer"]]),
                   absError = getOption("adestr_absError_outer", default = .adestr_options[["adestr_absError_outer"]]),
                   exact = FALSE,
                   early_futility_part = TRUE,
                   continuation_part = TRUE,
                   early_efficacy_part = TRUE,
                   conditional_integral = FALSE) {
  n1 <- n1(design, round = FALSE)
  if (exact){
    n1 <- ceiling(n1)
  }
  n2tmp <- n2_extrapol(design, seq(design@c1f, design@c1e, length.out=100))
  n2_min <- min(n2tmp)
  n2_max <- max(n2tmp)
  if (exact){
    n2_min <- ceiling(min(n2tmp))
    n2_max <- ceiling(max(n2tmp))
  }
  se1 <- sigma * sqrt((1L + two_armed) / n1)
  se2_n2min <- sigma * sqrt((1L + two_armed) / n2_min)
  se2_n2max <- sigma * sqrt((1L + two_armed) / n2_max)
  se2s <- c(se2_n2min, se2_n2max)
  infq <- if (absError==0) min(1e-6, tol) else min(1e-6, absError)
  minusinf_norm_1 <- get_norm_inf(infq, means = mu/se1, left = TRUE)
  plusinf_norm_1 <- get_norm_inf(infq, means = mu/se1, left = FALSE)
  minusinf_norm_2 <- get_norm_inf(infq, means = mu/se2s, left = TRUE)
  plusinf_norm_2 <- get_norm_inf(infq, means = mu/se2s, left = FALSE)

  futility_integral <- continuation_integral <- efficacy_integral <- NULL
  denom <- 1L
  if (conditional_integral) {
    denom <- 0
    if (early_futility_part)
      denom <- denom + pnorm(design@c1f, mean = mu/se1)
    if (continuation_part)
      denom <- denom + pnorm(design@c1e, mean = mu/se1) - pnorm(design@c1f, mean = mu/se1)
    if (early_efficacy_part)
      denom <- denom + 1 - pnorm(design@c1e, mean = mu/se1)
  }
  if (early_futility_part || early_efficacy_part) {
    mint1 <- \(x) {
      n2 <- n2_extrapol(design, x[1L, ])
      g1val <- g1(
        design = design,
        smean1 = x[1L, ] * se1,
        n1 = n1,
        mu = mu,
        sigma = sigma,
        two_armed = two_armed
      )
      sign(g1val) * exp(mlogf1_kv(x, n1, mu, sigma, two_armed) +
                          log(abs(g1val)))
    }
    if (early_futility_part) {
      futility_integral <- .hcubature(
        f = mint1,
        lowerLimit = min(minusinf_norm_1, design@c1f),
        upperLimit = design@c1f,
        tol = tol,
        maxEval = maxEval,
        absError = absError,
        vectorInterface = TRUE)
    }
    if (early_efficacy_part) {
      efficacy_integral <- .hcubature(
        f = mint1,
        lowerLimit = design@c1e,
        upperLimit = max(plusinf_norm_1, design@c1e),
        tol = tol,
        maxEval = maxEval,
        absError = absError,
        vectorInterface = TRUE)
    }
  }
  if (continuation_part) {
    if (!exact){
      mint2 <- \(x) {
        n2 <- n2_extrapol(design, x[1L, ])
        g2val <- g2(
          design = design,
          smean1 = x[1L, ] * se1,
          smean2 = x[2L, ] * sigma * sqrt((1L + two_armed) / n2),
          n1 = n1,
          n2 = n2,
          mu = mu,
          sigma = sigma,
          two_armed = two_armed)
        sign(g2val) * exp(mlogf2_kv(x, n1, n2, mu, sigma, two_armed) +
                            log(abs(g2val)))
      }
      continuation_integral <- .hcubature(
        f = mint2,
        lowerLimit = c(design@c1f, minusinf_norm_2),
        upperLimit = c(design@c1e, plusinf_norm_2),
        tol = tol,
        maxEval = maxEval,
        absError = absError,
        vectorInterface = TRUE)
    } else {
      mint2_exact <- \(x, n2, preimage) {
        g2val <- g2(
          design = design,
          smean1 = x[1L, ] * se1,
          smean2 = x[2L, ] * sigma * sqrt((1L + two_armed) / n2),
          n1 = n1,
          n2 = n2,
          mu = mu,
          sigma = sigma,
          two_armed = two_armed,
          sample_space = sample_space,
          preimage = preimage)
        sign(g2val) * exp(mlogf2_kv(x, n1, n2, mu, sigma, two_armed) +
                            log(abs(g2val)))
      }
      sample_space <- n2_preimage(design = design, sigma = sigma, two_armed = two_armed, smean_scale = FALSE)
      intexact <- function(region) {
        .hcubature(
          f = mint2_exact,
          lowerLimit = c(region$preimage[1L], minusinf_norm_2),
          upperLimit = c(region$preimage[2L], plusinf_norm_2),
          tol = tol,
          maxEval = maxEval,
          absError = absError,
          vectorInterface = TRUE,
          n2 = region$n2,
          preimage = region$preimage
        )
      }
      continuation_integral_list <- lapply(sample_space, intexact)
      continuation_integral <-
        as.list(sapply(X = names(continuation_integral_list[[1L]]), USE.NAMES = TRUE,
                       FUN = \(x) sum(sapply(continuation_integral_list, \(y)y[[x]]))))
    }
  }
  if (!early_efficacy_part && !early_futility_part && !continuation_part)
    return(list(overall_integral = list(integral = 0)))
  else {
    integrals <- list(futility_integral, continuation_integral, efficacy_integral)
    integrals <- integrals[lengths(integrals)!=0L]
    overall_integral <- as.list(sapply(X = names(integrals[[1L]]), USE.NAMES = TRUE,
                                       FUN = \(x) sum(sapply(integrals, \(y)y[[x]]))))
    overall_integral$integral <- overall_integral$integral / denom
    return(
      list(
        futility_integral = futility_integral,
        continuation_integral = continuation_integral,
        efficacy_integral = efficacy_integral,
        overall_integral = overall_integral
      )
    )
  }
}



# integrate distribution of statistic for 2 unknown parameter (mu, sigma^2) over adoptr design
#' @importFrom cubature hcubature
int_uv <- function(design,
                   g1 = \(...) 1L,
                   g2 = \(...) 1L,
                   mu,
                   sigma,
                   two_armed,
                   tol = getOption("adestr_tol_outer", default = .adestr_options[["adestr_tol_outer"]]),
                   maxEval = getOption("adestr_maxEval_outer", default = .adestr_options[["adestr_maxEval_outer"]]),
                   absError = getOption("adestr_absError_outer", default = .adestr_options[["adestr_absError_outer"]]),
                   exact = FALSE,
                   early_futility_part = TRUE,
                   continuation_part = TRUE,
                   early_efficacy_part = TRUE,
                   conditional_integral = FALSE) {
  n1 <- n1(design, round = FALSE)
  if (exact){
    n1 <- ceiling(n1)
  }
  se1 <- sigma * sqrt((1L + two_armed) / n1)
  v <- sigma^2
  df1 <- (1L + two_armed) * (n1-1L)
  n2tmp <- n2_extrapol(design, seq(design@c1f, design@c1e, length.out=100))
  n2_min <- min(n2tmp)
  n2_max <- max(n2tmp)
  if (exact){
    n2_min <- ceiling(min(n2tmp))
    n2_max <- ceiling(max(n2tmp))
  }
  df2_min <- (1L + two_armed) * (n2_min-1L)
  df2_max <- (1L + two_armed) * (n2_max-1L)
  df2s <- c(df2_min, df2_max)
  se2_n2min <- sigma * sqrt((1L + two_armed) / n2_min)
  se2_n2max <- sigma * sqrt((1L + two_armed) / n2_max)
  se2s <- c(se2_n2min, se2_n2max)
  infq <- if (absError==0) min(1e-6, tol) else min(1e-6, absError)
  minusinf_t_1 <- get_t_inf(infq, dfs = df1, means = mu/se1, left = TRUE)
  plusinf_t_1 <- get_t_inf(infq, dfs = df1, means = mu/se1, left = FALSE)
  minusinf_t_2 <- get_t_inf(infq, dfs = df2s, means = mu/se2s, left = TRUE)
  plusinf_t_2 <- get_t_inf(infq, dfs = df2s, means = mu/se2s, left = FALSE)
  minusinf_chi_1 <- get_chi_inf(infq, dfs = df1, variance = v, left = TRUE)
  plusinf_chi_1 <- get_chi_inf(infq, dfs = df1, variance = v, left = FALSE)
  minusinf_chi_2 <- get_chi_inf(infq, dfs = df2s, variance = v, left = TRUE)
  plusinf_chi_2 <- get_chi_inf(infq, dfs = df2s, variance = v, left = FALSE)

  futility_integral <- continuation_integral <- efficacy_integral <- NULL
  denom <- 1L
  if (conditional_integral) {
    denom <- 0
    if (early_futility_part)
      denom <- denom + pt(design@c1f, df1, mu/se1)
    if (continuation_part)
      denom <- denom + pt(design@c1e, df1, mu/se1) - pt(design@c1f, df1, mu/se1)
    if (early_efficacy_part)
      denom <- denom + 1 - pt(design@c1e, df1, mu/se1)
  }
  if (early_futility_part || early_efficacy_part) {
    mint1 <- \(x) {
      g1val <- g1(
        design = design,
        smean1 = x[1L, ] * sqrt((x[2L, ] * (1L + two_armed)) / n1),
        svar1 = x[2L, ],
        n1 = n1,
        mu = mu,
        sigma = sigma,
        two_armed = two_armed
      )
      sign(g1val) * exp(mlogf1_uv(x, n1, mu, sigma, two_armed) +
                          log(abs(g1val)))
    }
    if (early_futility_part) {
      futility_integral <- .hcubature(
        f = mint1,
        lowerLimit = c(min(minusinf_t_1, design@c1f), minusinf_chi_1),
        upperLimit = c(design@c1f, plusinf_chi_1),
        tol = tol,
        maxEval = maxEval,
        absError = absError,
        vectorInterface = TRUE)
    }
    if (early_efficacy_part) {
      efficacy_integral <- .hcubature(
        f = mint1,
        lowerLimit = c(design@c1e, minusinf_chi_1),
        upperLimit = c(max(plusinf_t_1, design@c1e), plusinf_chi_1),
        tol = tol,
        maxEval = maxEval,
        absError = absError,
        vectorInterface = TRUE)
    }
  }
  if (continuation_part) {
    if (!exact){
      mint2 <- \(x) {
        n2 <- n2_extrapol(design, x[1L, ])
        g2val <- g2(
          design = design,
          smean1 = x[1L, ] * sqrt((x[2L, ] * (1L + two_armed)) / n1),
          svar1  = x[2L,],
          smean2 = x[3L, ] * sqrt((x[4L, ] * (1L + two_armed)) / n2),
          svar2  = x[4L,],
          n1 = n1,
          n2 = n2,
          mu = mu,
          sigma = sigma,
          two_armed = two_armed)
        sign(g2val) * exp(mlogf2_uv(x, n1, n2, mu, sigma, two_armed) +
                            log(abs(g2val)))
      }
      continuation_integral <- .hcubature(
        f = mint2,
        lowerLimit = c(design@c1f, minusinf_chi_1, minusinf_t_2, minusinf_chi_2),
        upperLimit = c(design@c1e, plusinf_chi_1, plusinf_t_2, plusinf_chi_2),
        tol = tol,
        maxEval = maxEval,
        absError = absError,
        vectorInterface = TRUE)
    } else {
      mint2_exact <- \(x, n2, preimage) {
        g2val <- g2(
          design = design,
          smean1 = x[1L, ] * sqrt((x[2L, ] * (1L + two_armed)) / n1),
          svar1  = x[2L,],
          smean2 = x[3L, ] * sqrt((x[4L, ] * (1L + two_armed)) / n2),
          svar2  = x[4L,],
          n1 = n1,
          n2 = n2,
          mu = mu,
          sigma = sigma,
          two_armed = two_armed,
          sample_space = sample_space,
          preimage = preimage)
        sign(g2val) * exp(mlogf2_uv(x, n1, n2, mu, sigma, two_armed) +
                            log(abs(g2val)))
      }
      sample_space <- n2_preimage(design = design, sigma = sigma, two_armed = two_armed, smean_scale = FALSE)
      intexact <- function(region) {
        .hcubature(
          f = mint2_exact,
          lowerLimit = c(region$preimage[1L], minusinf_chi_1, minusinf_t_2, minusinf_chi_2),
          upperLimit = c(region$preimage[2L], plusinf_chi_1, plusinf_t_2, plusinf_chi_2),
          tol = tol,
          maxEval = maxEval,
          absError = absError,
          vectorInterface = TRUE,
          n2 = region$n2,
          preimage = region$preimage
        )
      }
      continuation_integral_list <- lapply(sample_space, intexact)
      continuation_integral <-
        as.list(sapply(X = names(continuation_integral_list[[1L]]), USE.NAMES = TRUE,
                       FUN = \(x) sum(sapply(continuation_integral_list, \(y)y[[x]]))))
    }
  }
  if (!early_efficacy_part && !early_futility_part && !continuation_part)
    return(list(overall_integral = list(integral = 0)))
  else {
    integrals <- list(futility_integral, continuation_integral, efficacy_integral)
    integrals <- integrals[lengths(integrals)!=0L]
    overall_integral <- as.list(sapply(X = names(integrals[[1L]]), USE.NAMES = TRUE,
                                       FUN = \(x) sum(sapply(integrals, \(y)y[[x]]))))
    overall_integral$integral <- overall_integral$integral / denom
    return(
      list(
        futility_integral = futility_integral,
        continuation_integral = continuation_integral,
        efficacy_integral = efficacy_integral,
        overall_integral = overall_integral
      )
    )
  }
}

#' @importFrom cubature hcubature
int_kv_full <- function(design,
                        g1 = \(...) 1L,
                        g2 = \(...) 1L,
                        mu,
                        sigma,
                        tol = getOption("adestr_tol_outer", default = .adestr_options[["adestr_tol_outer"]]),
                        maxEval = getOption("adestr_maxEval_outer", default = .adestr_options[["adestr_maxEval_outer"]]),
                        absError = getOption("adestr_absError_outer", default = .adestr_options[["adestr_absError_outer"]]),
                        exact = FALSE,
                        early_futility_part = TRUE,
                        continuation_part = TRUE,
                        early_efficacy_part = TRUE,
                        conditional_integral = FALSE) {
  n1 <- n1(design, round = FALSE)
  if (exact){
    n1 <- ceiling(n1)
  }
  n2tmp <- n2_extrapol(design, seq(design@c1f, design@c1e, length.out=100))
  n2_min <- min(n2tmp)
  n2_max <- max(n2tmp)
  if (exact){
    n2_min <- ceiling(min(n2tmp))
    n2_max <- ceiling(max(n2tmp))
  }
  se1_onearm <- sigma / sqrt(n1)
  se1_twoarm <- sqrt(2) * se1_onearm
  se2_n2min_onearm <- sigma / sqrt(n2_min)
  se2_n2max_onearm <- sigma / sqrt(n2_max)
  se2_n2min_twoarm <- sqrt(2) * se2_n2min_onearm
  se2_n2max_twoarm <- sqrt(2) * se2_n2max_onearm
  se2s_onearm <- c(se2_n2min_onearm, se2_n2max_onearm)
  se2s_twoarm <- c(se2_n2min_twoarm, se2_n2max_twoarm)

  infq <- if (absError==0) min(1e-6, tol) else min(1e-6, absError)

  minusinf_norm_onearm <- get_norm_inf(infq, means = 0, left = TRUE)
  plusinf_norm_onearm <- get_norm_inf(infq, means = 0, left = FALSE)
  minusinf_norm_1_twoarm <- get_norm_inf(infq, means = mu/se1_twoarm, left = TRUE)
  plusinf_norm_1_twoarm <- get_norm_inf(infq, means = mu/se1_twoarm, left = FALSE)
  minusinf_norm_2_twoarm <- get_norm_inf(infq, means = mu/se2s_twoarm, left = TRUE)
  plusinf_norm_2_twoarm <- get_norm_inf(infq, means = mu/se2s_twoarm, left = FALSE)

  futility_integral <- continuation_integral <- efficacy_integral <- NULL
  denom <- 1L
  if (conditional_integral) {
    denom <- 0
    if (early_futility_part)
      denom <- denom + pnorm(design@c1f, mean = mu/se1_twoarm)
    if (continuation_part)
      denom <- denom + pnorm(design@c1e, mean = mu/se1_twoarm) - pnorm(design@c1f, mean = mu/se1_twoarm)
    if (early_efficacy_part)
      denom <- denom + 1 - pnorm(design@c1e, mean = mu/se1_twoarm)
  }
  if (early_futility_part || early_efficacy_part) {
    mint1 <- \(x) {
      n2 <- n2_extrapol(design, x[1L, ])
      smean1 <- x[1L, ] * se1_twoarm
      smean1C <- x[2L, ] * se1_onearm
      g1val <- g1(
        design = design,
        smean1 = smean1,
        smean1T = smean1 + smean1C,
        n1 = n1,
        mu = mu,
        sigma = sigma
      )
      sign(g1val) * exp(mlogf1_kv_full(x, n1, mu, sigma) +
                          log(abs(g1val)))
    }
    if (early_futility_part) {
      futility_integral <- .hcubature(
        f = mint1,
        lowerLimit = c(min(minusinf_norm_1_twoarm, design@c1f), minusinf_norm_onearm),
        upperLimit = c(design@c1f, plusinf_norm_onearm),
        tol = tol,
        maxEval = maxEval,
        absError = absError,
        vectorInterface = TRUE)
    }
    if (early_efficacy_part) {
      efficacy_integral <- .hcubature(
        f = mint1,
        lowerLimit = c(design@c1e, minusinf_norm_onearm),
        upperLimit = c(max(plusinf_norm_1_twoarm, design@c1e), plusinf_norm_onearm),
        tol = tol,
        maxEval = maxEval,
        absError = absError,
        vectorInterface = TRUE)
    }
  }
  if (continuation_part) {
    if (!exact){
      mint2 <- \(x) {
        n2 <- n2_extrapol(design, x[1L, ])
        smean1 <- x[1L, ] * se1_twoarm
        smean1C <- x[2L, ] * se1_onearm
        se2_onearm <- sigma / sqrt(n2)
        smean2 <- x[3L, ] * se2_onearm * sqrt(2)
        smean2C <- x[4L, ] * se2_onearm
        g2val <- g2(
          design = design,
          smean1 = smean1,
          smean1T = smean1 + smean1C,
          smean2 = smean2,
          smean2T = smean2 + smean2C,
          n1 = n1,
          n2 = n2,
          mu = mu,
          sigma = sigma)
        sign(g2val) * exp(mlogf2_kv_full(x, n1, n2, mu, sigma) +
                            log(abs(g2val)))
      }
      continuation_integral <- .hcubature(
        f = mint2,
        lowerLimit = c(design@c1f, minusinf_norm_onearm, minusinf_norm_2_twoarm, minusinf_norm_onearm),
        upperLimit = c(design@c1e, plusinf_norm_onearm, plusinf_norm_2_twoarm, plusinf_norm_onearm),
        tol = tol,
        maxEval = maxEval,
        absError = absError,
        vectorInterface = TRUE)
    } else {
      mint2_exact <- \(x, n2, preimage) {
        smean1 <- x[1L, ] * se1_twoarm
        smean1C <- x[2L, ] * se1_onearm
        se2_onearm <- sigma / sqrt(n2)
        smean2 <- x[3L, ] * se2_onearm * sqrt(2)
        smean2C <- x[4L, ] * se2_onearm
        g2val <- g2(
          design = design,
          smean1 = smean1,
          smean1T = smean1 + smean1C,
          smean2 = smean2,
          smean2T = smean2 + smean2C,
          n1 = n1,
          n2 = n2,
          mu = mu,
          sigma = sigma,
          sample_space = sample_space,
          preimage = preimage)
        sign(g2val) * exp(mlogf2_kv_full(x, n1, n2, mu, sigma) +
                            log(abs(g2val)))
      }
      sample_space <- n2_preimage(design = design, sigma = sigma, two_armed = TRUE, smean_scale = FALSE)
      intexact <- function(region) {
        .hcubature(
          f = mint2_exact,
          lowerLimit = c(region$preimage[1L], minusinf_norm_onearm, minusinf_norm_2_twoarm, minusinf_norm_onearm),
          upperLimit = c(region$preimage[2L], plusinf_norm_onearm, plusinf_norm_2_twoarm, plusinf_norm_onearm),
          tol = tol,
          maxEval = maxEval,
          absError = absError,
          vectorInterface = TRUE,
          n2 = region$n2,
          preimage = region$preimage
        )
      }
      continuation_integral_list <- lapply(sample_space, intexact)
      continuation_integral <-
        as.list(sapply(X = names(continuation_integral_list[[1L]]), USE.NAMES = TRUE,
                       FUN = \(x) sum(sapply(continuation_integral_list, \(y)y[[x]]))))
    }
  }
  if (!early_efficacy_part && !early_futility_part && !continuation_part)
    return(list(overall_integral = list(integral = 0)))
  else {
    integrals <- list(futility_integral, continuation_integral, efficacy_integral)
    integrals <- integrals[lengths(integrals)!=0L]
    overall_integral <- as.list(sapply(X = names(integrals[[1L]]), USE.NAMES = TRUE,
                                       FUN = \(x) sum(sapply(integrals, \(y)y[[x]]))))
    overall_integral$integral <- overall_integral$integral / denom
    return(
      list(
        futility_integral = futility_integral,
        continuation_integral = continuation_integral,
        efficacy_integral = efficacy_integral,
        overall_integral = overall_integral
      )
    )
  }
}

#' @importFrom cubature hcubature
int_uv_full <- function(design,
                        g1 = \(...) 1L,
                        g2 = \(...) 1L,
                        mu,
                        sigma,
                        tol = getOption("adestr_tol_outer", default = .adestr_options[["adestr_tol_outer"]]),
                        maxEval = getOption("adestr_maxEval_outer", default = .adestr_options[["adestr_maxEval_outer"]]),
                        absError = getOption("adestr_absError_outer", default = .adestr_options[["adestr_absError_outer"]]),
                        exact = FALSE,
                        early_futility_part = TRUE,
                        continuation_part = TRUE,
                        early_efficacy_part = TRUE,
                        conditional_integral = FALSE) {
  n1 <- n1(design, round = FALSE)
  if (exact){
    n1 <- ceiling(n1)
  }
  se1_onearm <- sigma / sqrt(n1)
  se1_twoarm <- sqrt(2) * se1_onearm
  v <- sigma^2
  df1 <- 2L * (n1-1L)
  n2tmp <- n2_extrapol(design, seq(design@c1f, design@c1e, length.out=100))
  n2_min <- min(n2tmp)
  n2_max <- max(n2tmp)
  if (exact){
    n2_min <- ceiling(min(n2tmp))
    n2_max <- ceiling(max(n2tmp))
  }
  df2_min <- 2L * (n2_min-1L)
  df2_max <- 2L * (n2_max-1L)
  df2s <- c(df2_min, df2_max)
  se2_n2min_twoarm <- sigma * sqrt(2L / n2_min)
  se2_n2max_twoarm <- sigma * sqrt(2L / n2_max)

  se2s_twoarm <- c(se2_n2min_twoarm, se2_n2max_twoarm)
  infq <- if (absError==0) min(1e-6, tol) else min(1e-6, absError)

  minusinf_t_1_onearm <- get_t_inf(infq, dfs = df1, means = 0, left = TRUE)
  plusinf_t_1_onearm <- get_t_inf(infq, dfs = df1, means = 0, left = FALSE)
  minusinf_t_2_onearm <- get_t_inf(infq, dfs = df2s, means = 0, left = TRUE)
  plusinf_t_2_onearm <- get_t_inf(infq, dfs = df2s, means = 0, left = FALSE)

  minusinf_t_1_twoarm <- get_t_inf(infq, dfs = df1, means = mu/se1_twoarm, left = TRUE)
  plusinf_t_1_twoarm <- get_t_inf(infq, dfs = df1, means = mu/se1_twoarm, left = FALSE)
  minusinf_t_2_twoarm <- get_t_inf(infq, dfs = df2s, means = mu/se2s_twoarm, left = TRUE)
  plusinf_t_2_twoarm <- get_t_inf(infq, dfs = df2s, means = mu/se2s_twoarm, left = FALSE)

  minusinf_chi_1 <- get_chi_inf(infq, dfs = df1, variance = v, left = TRUE)
  plusinf_chi_1 <- get_chi_inf(infq, dfs = df1, variance = v, left = FALSE)
  minusinf_chi_2 <- get_chi_inf(infq, dfs = df2s, variance = v, left = TRUE)
  plusinf_chi_2 <- get_chi_inf(infq, dfs = df2s, variance = v, left = FALSE)


  futility_integral <- continuation_integral <- efficacy_integral <- NULL
  denom <- 1L
  if (conditional_integral) {
    denom <- 0
    if (early_futility_part)
      denom <- denom + pt(design@c1f, df1, mu/se1_twoarm)
    if (continuation_part)
      denom <- denom + pt(design@c1e, df1, mu/se1_twoarm) - pt(design@c1f, df1, mu/se1_twoarm)
    if (early_efficacy_part)
      denom <- denom + 1 - pt(design@c1e, df1, mu/se1_twoarm)
  }
  if (early_futility_part || early_efficacy_part) {
    mint1 <- \(x) {
      sdest1_onearm <- sqrt((x[3L, ]) / n1)
      smean1 <- x[1L, ] * sdest1_onearm * .SQRT2
      smean1C <- x[2L, ] * sdest1_onearm
      g1val <- g1(
        design = design,
        smean1 = smean1,
        smean1T = smean1 + smean1C,
        svar1 = x[3L, ],
        n1 = n1,
        mu = mu,
        sigma = sigma)
      sign(g1val) * exp(mlogf1_uv_full(x, n1, mu, sigma) +
                          log(abs(g1val)))
    }
    if (early_futility_part) {
      futility_integral <- .hcubature(
        f = mint1,
        lowerLimit = c(min(minusinf_t_1_twoarm, design@c1f), minusinf_t_1_onearm, minusinf_chi_1),
        upperLimit = c(design@c1f, plusinf_t_1_onearm, plusinf_chi_1),
        tol = tol,
        maxEval = maxEval,
        absError = absError,
        vectorInterface = TRUE)
    }
    if (early_efficacy_part) {
      efficacy_integral <- .hcubature(
        f = mint1,
        lowerLimit = c(design@c1e, minusinf_t_1_onearm, minusinf_chi_1),
        upperLimit = c(max(plusinf_t_1_twoarm, design@c1e), plusinf_t_1_onearm, plusinf_chi_1),
        tol = tol,
        maxEval = maxEval,
        absError = absError,
        vectorInterface = TRUE)
    }
  }
  if (continuation_part) {
    if (!exact){
      mint2 <- \(x) {
        n2 <- n2_extrapol(design, x[1L, ])
        sdest1_onearm <- sqrt(x[3L,] / n1)
        sdest2_onearm <- sqrt(x[6L,] / n2)
        smean1 <- x[1L, ] * sdest1_onearm * .SQRT2
        smean1C <- x[2L, ] * sdest1_onearm
        smean2 <- x[4L, ] * sdest2_onearm * .SQRT2
        smean2C <- x[5L, ] * sdest2_onearm
        g2val <- g2(
          design = design,
          smean1 = smean1,
          smean1T = smean1 + smean1C,
          svar1  = x[3L,],
          smean2 = smean2,
          smean2T = smean2 + smean2C,
          svar2  = x[6L,],
          n1 = n1,
          n2 = n2,
          mu = mu,
          sigma = sigma)
        sign(g2val) * exp(mlogf2_uv_full(x, n1, n2, mu, sigma) +
                            log(abs(g2val)))
      }
      continuation_integral <- .hcubature(
        f = mint2,
        lowerLimit = c(design@c1f, minusinf_t_1_onearm, minusinf_chi_1, minusinf_t_2_twoarm, minusinf_t_2_onearm , minusinf_chi_2),
        upperLimit = c(design@c1e, plusinf_t_1_onearm,  plusinf_chi_1,  plusinf_t_2_twoarm, plusinf_t_2_onearm,  plusinf_chi_2),
        tol = tol,
        maxEval = maxEval,
        absError = absError,
        vectorInterface = TRUE)
    } else {
      mint2_exact <- \(x, n2, preimage) {
        sdest1_onearm <- sqrt(x[3L,] / n1)
        sdest2_onearm <- sqrt(x[6L,] / n2)
        smean1 <- x[1L, ] * sdest1_onearm * .SQRT2
        smean1C <- x[2L, ] * sdest1_onearm
        smean2 <- x[4L, ] * sdest2_onearm * .SQRT2
        smean2C <- x[5L, ] * sdest2_onearm
        g2val <- g2(
          design = design,
          smean1 = smean1,
          smean1T = smean1 + smean1C,
          svar1  = x[3L,],
          smean2 = smean2,
          smean2T = smean2 + smean2C,
          svar2  = x[6L,],
          n1 = n1,
          n2 = n2,
          mu = mu,
          sigma = sigma,
          sample_space = sample_space,
          preimage = preimage)
        sign(g2val) * exp(mlogf2_uv_full(x, n1, n2, mu, sigma) +
                            log(abs(g2val)))
      }
      sample_space <- n2_preimage(design = design, sigma = sigma, two_armed = TRUE, smean_scale = FALSE)
      intexact <- function(region) {
        .hcubature(
          f = mint2_exact,
          lowerLimit = c(region$preimage[1L], minusinf_t_1_onearm, minusinf_chi_1, minusinf_t_2_twoarm, minusinf_t_2_onearm , minusinf_chi_2),
          upperLimit = c(region$preimage[2L], plusinf_t_1_onearm,  plusinf_chi_1,  plusinf_t_2_twoarm, plusinf_t_2_onearm,  plusinf_chi_2),
          tol = tol,
          maxEval = maxEval,
          absError = absError,
          vectorInterface = TRUE,
          n2 = region$n2,
          preimage = region$preimage
        )
      }
      continuation_integral_list <- lapply(sample_space, intexact)
      continuation_integral <-
        as.list(sapply(X = names(continuation_integral_list[[1L]]), USE.NAMES = TRUE,
                       FUN = \(x) sum(sapply(continuation_integral_list, \(y)y[[x]]))))
    }
  }
  if (!early_efficacy_part && !early_futility_part && !continuation_part)
    return(list(overall_integral = list(integral = 0)))
  else {
    integrals <- list(futility_integral, continuation_integral, efficacy_integral)
    integrals <- integrals[lengths(integrals)!=0L]
    overall_integral <- as.list(sapply(X = names(integrals[[1L]]), USE.NAMES = TRUE,
                                       FUN = \(x) sum(sapply(integrals, \(y)y[[x]]))))
    overall_integral$integral <- overall_integral$integral / denom
    return(
      list(
        futility_integral = futility_integral,
        continuation_integral = continuation_integral,
        efficacy_integral = efficacy_integral,
        overall_integral = overall_integral
      )
    )
  }
}

