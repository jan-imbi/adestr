setGeneric("dsmean", \(data_distribution, design, smean, mu, sigma, combine_components = TRUE,
                       tol = getOption("adestr_tol_outer", default = .adestr_options[["adestr_tol_outer"]]),
                       maxEval = getOption("adestr_maxEval_outer", default = .adestr_options[["adestr_maxEval_outer"]]),
                       absError = getOption("adestr_absError_outer", default = .adestr_options[["adestr_absError_outer"]]),
                       exact = FALSE, ...) standardGeneric("dsmean"))

#' @importFrom cubature hcubature
setMethod("dsmean", signature("Normal"),
          function(data_distribution, design, smean, mu, sigma, combine_components, tol, maxEval, absError, exact, ...){
            two_armed <- data_distribution@two_armed
            n1 <- n1(design, round = FALSE)
            if (exact){
              n1 <- ceiling(n1)
            }
            se1 <- sigma_to_se(sigma, n1, two_armed)
            zf <- design@c1f
            ze <- design@c1e
            if (exact) {
              sample_space <- n2_preimage(design = design, sigma = sigma, two_armed = two_armed, smean_scale = FALSE)
              ret <- sapply(smean,
                            \(x) {
                              partials <- sapply(
                                sample_space, \(region) {
                                  .hcubature(
                                    \(y) {
                                      n2 <- region$n2
                                      smean1 <- z_to_smean(y, n1, sigma, two_armed)
                                      se2 <- sigma_to_se(sigma, n2, two_armed)
                                      n <- (n1 + n2)
                                      n2se2 <- n2*se2
                                      z2 <- (n * x - n1 * smean1) / n2se2
                                      mf2_kv(rbind(y, z2), n1, n2, mu, sigma, two_armed) * (n / n2se2)
                                    },
                                    lowerLimit = region$preimage[[1L]],
                                    upperLimit = region$preimage[[2L]],
                                    tol = tol,
                                    absError = absError,
                                    maxEval = maxEval,
                                    vectorInterface = TRUE)$integral
                                })
                              z1 <- x / se1
                              futility_mass <- if (z1 < zf) f1_kv(z1, n1, mu, sigma, two_armed) / se1 else 0
                              efficacy_mass <- if (z1 > ze) f1_kv(z1, n1, mu, sigma, two_armed) / se1 else 0
                              if (combine_components)
                                return(sum(partials) + futility_mass + efficacy_mass)
                              else
                                return(c(partials,
                                         "efficacy" = efficacy_mass,
                                         "futility" = futility_mass))
                            }
              )
            } else {
              ret <-
                sapply(smean,
                       \(x) {
                         continuation <-
                           .hcubature(
                             \(y) {
                               n2 <- n2_extrapol(design, y)
                               smean1 <- z_to_smean(y, n1, sigma, two_armed)
                               se2 <- sigma_to_se(sigma, n2, two_armed)
                               n <- (n1 + n2)
                               n2se2 <- n2*se2
                               z2 <- (n * x - n1 * smean1) / n2se2
                               mf2_kv(rbind(y, z2), n1, n2, mu, sigma, two_armed) * (n / n2se2)
                             },
                             lowerLimit = zf,
                             upperLimit = ze,
                             tol = tol,
                             absError = absError,
                             maxEval = maxEval,
                             vectorInterface = TRUE)$integral
                         z1 <- x / se1
                         futility_mass <- if (z1 < zf) f1_kv(z1, n1, mu, sigma, two_armed) / se1 else 0
                         efficacy_mass <- if (z1 > ze) f1_kv(z1, n1, mu, sigma, two_armed) / se1 else 0
                         if (combine_components)
                           return(sum(continuation) + futility_mass + efficacy_mass)
                         else
                           return(c("continuation" = continuation,
                                    "efficacy" = efficacy_mass,
                                    "futility" = futility_mass))
                       })
            }
            if (combine_components)
              return(ret)
            else
              return(asplit(ret, MARGIN = 1))
          })

#' @importFrom cubature hcubature
setMethod("dsmean", signature("Student"),
          function(data_distribution, design, smean, mu, sigma, combine_components, tol, maxEval, absError, exact, ...){
            two_armed <- data_distribution@two_armed
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
            if (exact) {
              sample_space <- n2_preimage(design = design, sigma = sigma, two_armed = two_armed, smean_scale = FALSE)
              ret <- sapply(smean,
                            \(x) {
                              partials <- sapply(
                                sample_space, \(region) {
                                  .hcubature(
                                    \(y) {
                                      n2 <- region$n2
                                      smean1 <- t_to_smean(y[1L,,drop=FALSE], y[2L,,drop=FALSE], n1, two_armed)
                                      se2hat <- sqrt(y[3L,,drop=FALSE] * (1L + two_armed) / n2)
                                      n <- (n1 + n2)
                                      n2se2hat <- n2 * se2hat
                                      t2 <- (n * x - n1 * smean1) / n2se2hat
                                      mf2_uv(rbind(y[1L:2L,,drop=FALSE], t2, y[3L,,drop=FALSE]), n1, n2, mu, sigma, two_armed) * (n / n2se2hat)
                                    },
                                    lowerLimit = c(region$preimage[[1L]],
                                                   minusinf_chi_1,
                                                   minusinf_chi_2),
                                    upperLimit = c(region$preimage[[2L]],
                                                   plusinf_chi_1,
                                                   plusinf_chi_2),
                                    tol = tol,
                                    absError = absError,
                                    maxEval = maxEval,
                                    vectorInterface = TRUE)$integral
                                })
                              futility_mass <- .hcubature(
                                f = \(y) {
                                  z <- rbind(y[1L,,drop=FALSE], (x / (y[1L,,drop=FALSE] * (1L + two_armed)) )^2 * n1)
                                  mf1_uv(z, n1, mu, sigma, two_armed) * (sign(y[1L,])==sign(x)) * 2L *  (abs(x) / (y[1L,] * (1L + two_armed))^2) * n1
                                },
                                lowerLimit = if (x<0) min(design@c1f, minusinf_t_1) else 0,
                                upperLimit = if (x<0 & design@c1f > 0) 0 else design@c1f,
                                tol = tol,
                                absError = absError,
                                maxEval = maxEval,
                                vectorInterface = TRUE)$integral
                              efficacy_mass <- .hcubature(
                                f = \(y) {
                                  z <- rbind(y[1L,,drop=FALSE], (x / (y[1L,,drop=FALSE] * (1L + two_armed)) )^2 * n1)
                                  mf1_uv(z, n1, mu, sigma, two_armed) * (sign(y[1L,])==sign(x)) * 2L *  (abs(x) / (y[1L,] * (1L + two_armed))^2) * n1
                                },
                                lowerLimit = if (x > 0 & design@c1e < 0) 0 else design@c1e,
                                upperLimit = if (x > 0) max(design@c1e, plusinf_t_1) else 0,
                                tol = tol,
                                absError = absError,
                                maxEval = maxEval,
                                vectorInterface = TRUE)$integral
                              if (combine_components)
                                return(sum(partials) + futility_mass + efficacy_mass)
                              else
                                return(c(partials, "futility" =  futility_mass , "efficacy" =  efficacy_mass))
                            })
            } else {
              ret <- sapply(smean,
                            \(x) {
                              continuation <- .hcubature(
                                \(y) {
                                  n2 <- n2_extrapol(design, y[1L,,drop=FALSE])
                                  smean1 <- t_to_smean(y[1L,,drop=FALSE], y[2L,,drop=FALSE], n1, two_armed)
                                  se2hat <- sqrt(y[3L,,drop=FALSE] * (1L + two_armed) / n2)
                                  n <- (n1 + n2)
                                  n2se2hat <- n2 * se2hat
                                  t2 <- (n * x - n1 * smean1) / n2se2hat
                                  mf2_uv(rbind(y[1L:2L,,drop=FALSE], t2, y[3L,,drop=FALSE]), n1, n2, mu, sigma, two_armed) * (n / n2se2hat)
                                },
                                lowerLimit = c(design@c1f,
                                               minusinf_chi_1,
                                               minusinf_chi_2),
                                upperLimit = c(design@c1e,
                                               plusinf_chi_1,
                                               plusinf_chi_2),
                                tol = tol,
                                absError = absError,
                                maxEval = maxEval,
                                vectorInterface = TRUE)$integral
                              cond_dens <- \(y) {
                                z <- rbind(y[1L,,drop=FALSE], (x / (y[1L,,drop=FALSE] * (1L + two_armed)) )^2 * n1)
                                mf1_uv(z, n1, mu, sigma, two_armed) * (sign(y[1L,])==sign(x)) * 2L *  (abs(x) / (y[1L,] * (1L + two_armed))^2) * n1
                              }
                              ll_fut <- if (x<0) min(design@c1f, minusinf_t_1) else 0
                              ul_fut <- if (x<0 & design@c1f > 0) 0 else design@c1f
                              ll_eff <- if (x > 0 & design@c1e < 0) 0 else design@c1e
                              ul_eff <- if (x > 0) max(design@c1e, plusinf_t_1) else 0
                              testseq <- seq(ll_fut, ul_fut, length.out=15)
                              testseq <- testseq[testseq!=0]
                              if (!(abs(x-mu) * sqrt(n1 / (1L + two_armed))/sigma > qnorm(infq, lower.tail = FALSE))){
                                if (!any(cond_dens(t(testseq)) >1e-3)) { # TODO: maybe look for better magic number
                                  testseq <- seq(ll_fut, ul_fut, length.out=1e4)
                                  testseq <- testseq[testseq!=0]
                                  cd <- cond_dens(t(testseq))
                                  if (any(cd>1e-3)) {
                                    midpoint <- testseq[which(cd>1e-3)[1]]
                                    ll_fut <- ul_fut - 2* (ul_fut - midpoint)
                                  }
                                }
                              }
                              futility_mass <- .hcubature(
                                f = cond_dens,
                                lowerLimit = ll_fut,
                                upperLimit = ul_fut,
                                tol = tol,
                                absError = absError,
                                maxEval = maxEval,
                                vectorInterface = TRUE)$integral
                              efficacy_mass <- .hcubature(
                                f = cond_dens,
                                lowerLimit = if (x > 0 & design@c1e < 0) 0 else design@c1e,
                                upperLimit = if (x > 0) max(design@c1e, plusinf_t_1) else 0,
                                tol = tol,
                                absError = absError,
                                maxEval = maxEval,
                                vectorInterface = TRUE)$integral
                              # browser()
                              if (combine_components)
                                return(futility_mass + efficacy_mass + continuation)
                              else
                                return(c("continuation" = continuation, "futility" =  futility_mass , "efficacy" =  efficacy_mass))
                            }
              )
              if (combine_components)
                return(ret)
              else
                return(asplit(ret, MARGIN = 1))
            }
          })

setGeneric("dsmeanT", \(data_distribution, design, smeanT, mu, sigma, combine_components = TRUE,
                       tol = getOption("adestr_tol_outer", default = .adestr_options[["adestr_tol_outer"]]),
                       maxEval = getOption("adestr_maxEval_outer", default = .adestr_options[["adestr_maxEval_outer"]]),
                       absError = getOption("adestr_absError_outer", default = .adestr_options[["adestr_absError_outer"]]),
                       exact = FALSE, ...) standardGeneric("dsmeanT"))


setMethod("dsmeanT", signature("Normal"),
          function(data_distribution, design, smeanT, mu, sigma, combine_components = TRUE,
                    tol = getOption("adestr_tol_outer", default = .adestr_options[["adestr_tol_outer"]]),
                    maxEval = getOption("adestr_maxEval_outer", default = .adestr_options[["adestr_maxEval_outer"]]),
                    absError = getOption("adestr_absError_outer", default = .adestr_options[["adestr_absError_outer"]]),
                    exact = FALSE, ...) {
  if (!data_distribution@two_armed) {
    warning(paste0("dsmeanT calculates the density function of the overall sample mean for the treatment group in a two-armed trial. ",
                   "The specified data_distribution object is one-armed, which will be ignored. Use dsmean if you want to calculate the density for a one-armed trial."))
  }
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


  if (exact) {
    sample_space <- n2_preimage(design = design, sigma = sigma, two_armed = TRUE, smean_scale = FALSE)
    ret <- sapply(smeanT,
                  \(x) {
                    partials <- sapply(
                      sample_space, \(region) {
                        .hcubature(
                          \(y){
                            n2 <- region$n2
                            smean1 <- y[1L,,drop=FALSE] * se1_twoarm
                            smean1C <- y[2L,,drop=FALSE] * se1_onearm
                            se2_onearm <- sigma / sqrt(n2)
                            smean2 <- y[3L,,drop=FALSE] * se2_onearm * sqrt(2)
                            n <- (n1 + n2)
                            smean <- (n1 * smean1 + n2 * smean2) / n
                            smeanC <- x - smean
                            smean2C <- (n * smeanC - n1 * smean1C) / n2
                            mf2_kv_full(rbind(y, smean2C / se2_onearm), n1, n2, mu, sigma) * (n / (n2 * se2_onearm))
                          },
                          lowerLimit = c(region$preimage[[1L]],
                                         minusinf_norm_onearm,
                                         minusinf_norm_2_twoarm),
                          upperLimit = c(region$preimage[[2L]],
                                         plusinf_norm_onearm,
                                         plusinf_norm_2_twoarm),
                          tol = tol,
                          maxEval = maxEval,
                          absError = absError,
                          vectorInterface = TRUE)$integral
                      })
                      futility_mass <- .hcubature(\(y) {
                        smean1 <- y[1L,,drop=FALSE] * se1_twoarm
                        smean1C <- x - smean1
                        mf1_kv_full(rbind(y, smean1C / se1_onearm), n1, mu, sigma) / se1_onearm
                      },
                      lowerLimit = c(min(minusinf_norm_1_twoarm, design@c1f)),
                      upperLimit = c(design@c1f),
                      tol = tol,
                      maxEval = maxEval,
                      absError = absError,
                      vectorInterface = TRUE)$integral
                      efficacy_mass <- .hcubature(\(y) {
                        smean1 <- y[1L,,drop=FALSE] * se1_twoarm
                        smean1C <- x - smean1
                        mf1_kv_full(rbind(y, smean1C / se1_onearm), n1, mu, sigma) / se1_onearm
                      },
                      lowerLimit = c(design@c1e),
                      upperLimit = c(max(plusinf_norm_1_twoarm, design@c1e)),
                      tol = tol,
                      maxEval = maxEval,
                      absError = absError,
                      vectorInterface = TRUE)$integral
                    if (combine_components)
                      return(sum(partials) + futility_mass + efficacy_mass)
                    else
                      return(c(partials, "futility" =  futility_mass , "efficacy" =  efficacy_mass))
                  })

  } else {
    ret <- sapply(smeanT, \(x) {
      continuation <-  .hcubature(\(y){
        n2 <- n2_extrapol(design, y[1L,,drop=FALSE])
        smean1 <- y[1L,,drop=FALSE] * se1_twoarm
        smean1C <- y[2L,,drop=FALSE] * se1_onearm
        se2_onearm <- sigma / sqrt(n2)
        smean2 <- y[3L,,drop=FALSE] * se2_onearm * sqrt(2)
        n <- (n1 + n2)
        smean <- (n1 * smean1 + n2 * smean2) / n
        smeanC <- x - smean
        smean2C <- (n * smeanC - n1 * smean1C) / n2
        mf2_kv_full(rbind(y, smean2C / se2_onearm), n1, n2, mu, sigma) * (n / (n2 * se2_onearm))
      },
      lowerLimit = c(design@c1f,
                     minusinf_norm_onearm,
                     minusinf_norm_2_twoarm),
      upperLimit = c(design@c1e,
                     plusinf_norm_onearm,
                     plusinf_norm_2_twoarm),
      tol = tol,
      maxEval = maxEval,
      absError = absError,
      vectorInterface = TRUE)$integral

      futility_mass <- .hcubature(\(y) {
        smean1 <- y[1L,,drop=FALSE] * se1_twoarm
        smean1C <- x - smean1
        mf1_kv_full(rbind(y, smean1C / se1_onearm), n1, mu, sigma) / se1_onearm
        },
        lowerLimit = c(min(minusinf_norm_1_twoarm, design@c1f)),
        upperLimit = c(design@c1f),
        tol = tol,
        maxEval = maxEval,
        absError = absError,
        vectorInterface = TRUE)$integral
      efficacy_mass <- .hcubature(\(y) {
        smean1 <- y[1L,,drop=FALSE] * se1_twoarm
        smean1C <- x - smean1
        mf1_kv_full(rbind(y, smean1C / se1_onearm), n1, mu, sigma) / se1_onearm
      },
      lowerLimit = c(design@c1e),
      upperLimit = c(max(plusinf_norm_1_twoarm, design@c1e)),
      tol = tol,
      maxEval = maxEval,
      absError = absError,
      vectorInterface = TRUE)$integral
      if (combine_components)
        return(sum(continuation) + futility_mass + efficacy_mass)
      else
        return(c("continuation" = continuation, "futility" =  futility_mass , "efficacy" =  efficacy_mass))
    })
  }
  if (combine_components)
    return(ret)
  else
    return(asplit(ret, MARGIN = 1))
})

setMethod("dsmeanT", signature("Student"),
          function(data_distribution, design, smeanT, mu, sigma, combine_components = TRUE,
                    tol = getOption("adestr_tol_outer", default = .adestr_options[["adestr_tol_outer"]]),
                    maxEval = getOption("adestr_maxEval_outer", default = .adestr_options[["adestr_maxEval_outer"]]),
                    absError = getOption("adestr_absError_outer", default = .adestr_options[["adestr_absError_outer"]]),
                    exact = FALSE, ...) {
  if (!data_distribution@two_armed) {
    warning(paste0("dsmeanT calculates the density function of the overall sample mean for the treatment group in a two-armed trial. ",
            "The specified data_distribution object is one-armed, which will be ignored. Use dsmean if you want to calculate the density for a one-armed trial."))
  }
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

  if (exact) {
    sample_space <- n2_preimage(design = design, sigma = sigma, two_armed = TRUE, smean_scale = FALSE)
    ret <- sapply(smeanT,
                  \(x) {
                    partials <- sapply(
                      sample_space, \(region) {
                        .hcubature(
                          \(y){
                            n2 <- region$n2
                            sdest1_onearm <- sqrt(y[3L, ] / n1)
                            sdest2_onearm <- sqrt(y[5L, ] / n2)
                            smean1 <- y[1L, ] * sdest1_onearm * .SQRT2
                            smean1C <- y[2L, ] * sdest1_onearm
                            smean2 <- y[4L, ] * sdest2_onearm * .SQRT2
                            n <- n1 + n2
                            smean <- (n1 * smean1 + n2 * smean2) / n
                            smeanC <- x - smean
                            smean2C <- (n * smeanC - n1 * smean1C)/n2
                            mf2_uv_full(rbind(y[1L:4L,,drop=FALSE], smean2C / sdest2_onearm, y[5L,,drop=FALSE]), n1, n2, mu, sigma) * n / (n2 * sdest2_onearm)
                          },
                          lowerLimit = c(region$preimage[[1L]],
                                         minusinf_t_1_onearm,
                                         minusinf_chi_1,
                                         minusinf_t_2_twoarm,
                                         minusinf_chi_2
                          ),
                          upperLimit = c(region$preimage[[2L]],
                                         plusinf_t_1_onearm,
                                         plusinf_chi_1,
                                         plusinf_t_2_twoarm,
                                         plusinf_chi_2),
                          tol = tol,
                          maxEval = maxEval,
                          absError = absError,
                          vectorInterface = TRUE)$integral
                      })
                    futility_mass <- .hcubature(\(y) {
                      sdest1_onearm <- sqrt(y[2L,,drop=FALSE ] / n1)
                      smean1 <- y[1L,,drop=FALSE ] * sdest1_onearm * .SQRT2
                      smean1C <- x - smean1
                      mf1_uv_full(rbind(y[1L,,drop=FALSE], smean1C / sdest1_onearm, y[2L,,drop=FALSE]), n1, mu, sigma) / sdest1_onearm
                    },
                    lowerLimit = c(min(minusinf_t_1_twoarm, design@c1f), minusinf_chi_1),
                    upperLimit = c(design@c1f, plusinf_chi_1),
                    tol = tol,
                    maxEval = maxEval,
                    absError = absError,
                    vectorInterface = TRUE)$integral
                    efficacy_mass <- .hcubature(\(y) {
                      sdest1_onearm <- sqrt(y[2L,,drop=FALSE ] / n1)
                      smean1 <- y[1L,,drop=FALSE ] * sdest1_onearm * .SQRT2
                      smean1C <- x - smean1
                      mf1_uv_full(rbind(y[1L,,drop=FALSE], smean1C / sdest1_onearm, y[2L,,drop=FALSE]), n1, mu, sigma) / sdest1_onearm
                    },
                    lowerLimit = c(design@c1e, minusinf_chi_1),
                    upperLimit = c(max(plusinf_t_1_twoarm, design@c1e), plusinf_chi_1),
                    tol = tol,
                    maxEval = maxEval,
                    absError = absError,
                    vectorInterface = TRUE)$integral
                    if (combine_components)
                      return(sum(partials) + futility_mass + efficacy_mass)
                    else
                      return(c(partials, "futility" =  futility_mass , "efficacy" =  efficacy_mass))
                  })

  } else {
    ret <- sapply(smeanT, \(x) {
      continuation <-  .hcubature(\(y){
        n2 <- n2_extrapol(design, y[1L, ])
        sdest1_onearm <- sqrt(y[3L, ] / n1)
        sdest2_onearm <- sqrt(y[5L, ] / n2)
        smean1 <- y[1L, ] * sdest1_onearm * .SQRT2
        smean1C <- y[2L, ] * sdest1_onearm
        smean2 <- y[4L, ] * sdest2_onearm * .SQRT2
        n <- n1 + n2
        smean <- (n1 * smean1 + n2 * smean2) / n
        smeanC <- x - smean
        smean2C <- (n * smeanC - n1 * smean1C)/n2
        mf2_uv_full(rbind(y[1L:4L,,drop=FALSE], smean2C / sdest2_onearm, y[5L,,drop=FALSE]), n1, n2, mu, sigma) * n / (n2 * sdest2_onearm)
      },
      lowerLimit = c(design@c1f,
                     minusinf_t_1_onearm,
                     minusinf_chi_1,
                     minusinf_t_2_twoarm,
                     minusinf_chi_2
                     ),
      upperLimit = c(design@c1e,
                     plusinf_t_1_onearm,
                     plusinf_chi_1,
                     plusinf_t_2_twoarm,
                     plusinf_chi_2),
      tol = tol,
      maxEval = maxEval,
      absError = absError,
      vectorInterface = TRUE)$integral

      futility_mass <- .hcubature(\(y) {
        sdest1_onearm <- sqrt(y[2L,,drop=FALSE ] / n1)
        smean1 <- y[1L,,drop=FALSE ] * sdest1_onearm * .SQRT2
        smean1C <- x - smean1
        mf1_uv_full(rbind(y[1L,,drop=FALSE], smean1C / sdest1_onearm, y[2L,,drop=FALSE]), n1, mu, sigma) / sdest1_onearm
      },
      lowerLimit = c(min(minusinf_t_1_twoarm, design@c1f), minusinf_chi_1),
      upperLimit = c(design@c1f, plusinf_chi_1),
      tol = tol,
      maxEval = maxEval,
      absError = absError,
      vectorInterface = TRUE)$integral
      efficacy_mass <- .hcubature(\(y) {
        sdest1_onearm <- sqrt(y[2L,,drop=FALSE ] / n1)
        smean1 <- y[1L,,drop=FALSE ] * sdest1_onearm * .SQRT2
        smean1C <- x - smean1
        mf1_uv_full(rbind(y[1L,,drop=FALSE], smean1C / sdest1_onearm, y[2L,,drop=FALSE]), n1, mu, sigma) / sdest1_onearm
      },
      lowerLimit = c(design@c1e, minusinf_chi_1),
      upperLimit = c(max(plusinf_t_1_twoarm, design@c1e), plusinf_chi_1),
      tol = tol,
      maxEval = maxEval,
      absError = absError,
      vectorInterface = TRUE)$integral
      if (combine_components)
        return(sum(continuation) + futility_mass + efficacy_mass)
      else
        return(c("continuation" = continuation, "futility" =  futility_mass , "efficacy" =  efficacy_mass))
    })
  }
  if (combine_components)
    return(ret)
  else
    return(asplit(ret, MARGIN = 1))

})












