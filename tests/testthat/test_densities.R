library(cubature)
infq <- 1e-9
n1 <- 30
n2 <- 50
mu <- .3
sigma <- 2.1

v <- sigma^2
df1 <- n1-1L
df2 <- n2-1L
minusinf_norm_onearm <- qnorm(infq, mean = mu*sqrt(n1)/sigma, lower.tail = TRUE)
minusinf_norm_twoarm <- qnorm(infq, mean = mu*sqrt(n1/2)/sigma, lower.tail = TRUE)
minusinf_norm_control <- qnorm(infq, lower.tail = TRUE)
plusinf_norm_onearm <- qnorm(infq, mean = mu*sqrt(n1)/sigma, lower.tail = FALSE)
plusinf_norm_twoarm <- qnorm(infq, mean = mu*sqrt(n1/2)/sigma, lower.tail = FALSE)
plusinf_norm_control <- qnorm(infq, lower.tail = FALSE)

minusinf_chi_onearm_1 <- qchisq(infq, df1, lower.tail = TRUE) * v / df1
minusinf_chi_onearm_2 <- qchisq(infq, df2, lower.tail = TRUE) * v / df2
plusinf_chi_onearm_1  <- qchisq(infq, df1, lower.tail = FALSE) * v / df1
plusinf_chi_onearm_2  <- qchisq(infq, df2, lower.tail = FALSE) * v / df2

minusinf_chi_twoarm_1 <- qchisq(infq, 2*df1, lower.tail = TRUE) * v / (2L*df1)
minusinf_chi_twoarm_2 <- qchisq(infq, 2*df2, lower.tail = TRUE) * v / (2L*df2)
plusinf_chi_twoarm_1  <- qchisq(infq, 2*df1, lower.tail = FALSE) * v / (2L*df1)
plusinf_chi_twoarm_2  <- qchisq(infq, 2*df2, lower.tail = FALSE) * v / (2L*df2)

test_that("one-armed trial, known variance, integral of first-stage density equals 1",
          {
            expect_equal(
              hcubature(
                mf1_kv,
                lowerLimit = minusinf_norm_onearm,
                upperLimit = plusinf_norm_onearm,
                vectorInterface = TRUE,
                n1 = n1,
                mu = mu,
                sigma = sigma,
                two_armed = FALSE
              )$integral,
              1
            )
          })
test_that("one-armed trial, known variance, integral of second-stage density equals 1",
          {
            expect_equal(
              hcubature(
                mf2_kv,
                lowerLimit = rep(minusinf_norm_onearm, 2),
                upperLimit = rep(plusinf_norm_onearm, 2),
                vectorInterface = TRUE,
                n1 = n1,
                n2 = n2,
                mu = mu,
                sigma = sigma,
                two_armed = FALSE
              )$integral,
              1
            )
          })
test_that("two-armed trial, effect difference only, known variance, integral of first-stage density equals 1",
          {
            expect_equal(
              hcubature(
                mf1_kv,
                lowerLimit = minusinf_norm_twoarm,
                upperLimit = plusinf_norm_twoarm,
                vectorInterface = TRUE,
                n1 = n1,
                mu = mu,
                sigma = sigma,
                two_armed = TRUE)$integral,
              1
            )
          })
test_that("two-armed trial, effect difference only, known variance, integral of second-stage density equals 1",
          {
            expect_equal(
              hcubature(
                mf2_kv,
                lowerLimit = rep(minusinf_norm_twoarm, 2),
                upperLimit = rep(plusinf_norm_twoarm, 2),
                vectorInterface = TRUE,
                n1 = n1,
                n2 = n2,
                mu = mu,
                sigma = sigma,
                two_armed = TRUE)$integral,
              1,
              tolerance = 1e-6
            )
          })

test_that("two-armed trial, full sampling distribution, known variance, integral of first-stage density equals 1",
          {
            expect_equal(
              hcubature(
                mf1_kv_full,
                lowerLimit = c(minusinf_norm_onearm, minusinf_norm_twoarm),
                upperLimit = c(plusinf_norm_onearm, plusinf_norm_twoarm),
                vectorInterface = TRUE,
                n1 = n2,
                mu = mu,
                sigma = sigma)$integral,
              1,
              tolerance = 1e-5
            )
          })

test_that("two-armed trial, full sampling distribution, known variance, integral of second-stage density equals 1",
          {
            expect_equal(
              hcubature(
                mf2_kv_full,
                lowerLimit = rep(c(minusinf_norm_onearm, minusinf_norm_twoarm), 2),
                upperLimit = rep(c(plusinf_norm_onearm, plusinf_norm_twoarm), 2),
                vectorInterface = TRUE,
                n1 = n1,
                n2 = n2,
                mu = mu,
                sigma = sigma)$integral,
              1,
              tolerance = 1e-6
            )
          })

test_that("one-armed trial, unknown variance, integral of first-stage density equals 1",
          {
            expect_equal(
              hcubature(
                mf1_uv,
                lowerLimit = c(minusinf_norm_onearm, minusinf_chi_onearm_1),
                upperLimit = c(plusinf_norm_onearm, plusinf_chi_onearm_1),
                vectorInterface = TRUE,
                n1 = n1,
                mu = mu,
                sigma = sigma,
                two_armed = FALSE)$integral,
              1,
              tolerance = 1e-5
            )
          })

test_that("one-armed trial, unknown variance, integral of second-stage density equals 1",
          {
            expect_equal(
              hcubature(
                mf2_uv,
                lowerLimit = c(minusinf_norm_onearm, minusinf_chi_onearm_1, minusinf_norm_onearm, minusinf_chi_onearm_2),
                upperLimit = c(plusinf_norm_onearm, plusinf_chi_onearm_1, plusinf_norm_onearm, plusinf_chi_onearm_2),
                vectorInterface = TRUE,
                n1 = n1,
                n2 = n2,
                mu = mu,
                sigma = sigma,
                two_armed = FALSE)$integral,
              1,
              tolerance = 1e-5
            )
          })

test_that("two-armed trial, effect difference only, unknown variance, integral of first-stage density equals 1",
          {
            expect_equal(
              hcubature(
                mf1_uv,
                lowerLimit = c(minusinf_norm_twoarm, minusinf_chi_twoarm_1),
                upperLimit = c(plusinf_norm_twoarm, plusinf_chi_twoarm_1),
                vectorInterface = TRUE,
                n1 = n1,
                mu = mu,
                sigma = sigma,
                two_armed = TRUE)$integral,
              1,
              tolerance = 1e-6
            )
          })

test_that("two-armed trial, effect difference only, unknown variance, integral of second-stage density equals 1",
          {
            expect_equal(
              hcubature(
                mf2_uv,
                lowerLimit = c(minusinf_norm_twoarm, minusinf_chi_twoarm_1, minusinf_norm_twoarm, minusinf_chi_twoarm_2),
                upperLimit = c(plusinf_norm_twoarm, plusinf_chi_twoarm_1, plusinf_norm_twoarm, plusinf_chi_twoarm_2),
                vectorInterface = TRUE,
                n1 = n1,
                n2 = n2,
                mu = mu,
                sigma = sigma,
                two_armed = TRUE)$integral,
              1,
              tolerance = 1e-6
            )
          })

test_that("two-armed trial, full sampling distribution, unknown variance, integral of first-stage density equals 1",
          {
            expect_equal(
              hcubature(
                mf1_uv_full,
                lowerLimit = c(minusinf_norm_twoarm, minusinf_norm_control, minusinf_chi_twoarm_1),
                upperLimit = c(plusinf_norm_twoarm, plusinf_norm_control, plusinf_chi_twoarm_1),
                vectorInterface = TRUE,
                n1 = n1,
                mu = mu,
                sigma = sigma)$integral,
              1,
              tolerance = 1e-5
            )
          })

test_that("two-armed trial, full sampling distribution, unknown variance, integral of second-stage density equals 1",
          {
            expect_equal(
              hcubature(
                mf2_uv_full,
                lowerLimit = c(minusinf_norm_twoarm, minusinf_norm_control, minusinf_chi_twoarm_1, minusinf_norm_twoarm, minusinf_norm_control, minusinf_chi_twoarm_2),
                upperLimit = c(plusinf_norm_twoarm, plusinf_norm_control, plusinf_chi_twoarm_1, plusinf_norm_twoarm, plusinf_norm_control, plusinf_chi_twoarm_2),
                vectorInterface = TRUE,
                n1 = n1,
                n2 = n2,
                mu = mu,
                sigma = sigma,
                maxEval = 1e7)$integral,
              1,
              tolerance = 1e-3
            )
          })

