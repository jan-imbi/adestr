set.seed(321)
dat <- data.frame(
  endpoint = rnorm(sum(c(56, 56, 47, 47)), mean = rep(c(.3, 0, .3, 0), c(56, 56, 47, 47))),
  group = factor(rep(c("ctl", "trt", "ctl", "trt"), c(56,56,47,47))),
  stage = rep(c(1L, 2L), c(56*2, 47*2))
)
test_that("Analysis function doesn't throw an error.",
          {
            expect_error(analyze(
              data = dat,
              statistics = list(FirstStageSampleMean(), SampleMean(), NaiveCI()),
              data_distribution = Normal(TRUE),
              sigma = 1,
              design = get_example_design(TRUE)
            ),
            NA)
          })




