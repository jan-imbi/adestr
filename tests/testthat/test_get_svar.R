test_that("get_overall_svar_onearm correctly converts stage-wise sample means and variances to overall variance.",
          {
            n1 <- 10
            n2 <- 20
            x1 <- rnorm(n1)
            x2 <- rnorm(n2)

            smean1 <- mean(x1)
            smean2 <- mean(x2)
            svar1 <- var(x1)
            svar2 <- var(x2)

            expect_equal(get_overall_svar_onearm(smean1, svar1, smean2, svar2, n1, n2), var(c(x1,x2)))
          })

test_that("get_overall_svar_twoarm correctly converts stage-wise sample means and variances to overall variance.",
          {
            n1 <- 10
            n2 <- 20
            x1<- rnorm(n1)
            y1 <- rnorm(n1)
            x2<- rnorm(n2)
            y2 <- rnorm(n2)

            smean1T <- mean(x1)
            smean2T <- mean(x2)
            smean1C <- mean(y1)
            smean2C <- mean(y2)
            smean1 <- smean1T - smean1C
            smean2 <- smean2T - smean2C
            svar1 <- (var(x1) + var(y1))/2
            svar2 <- (var(x2) + var(y2))/2

            expect_equal(get_overall_svar_twoarm(smean1, smean1T, svar1, smean2, smean2T, svar2, n1, n2),
                         (var(c(x1, x2)) + var(c(y1, y2)))/2)
          })
