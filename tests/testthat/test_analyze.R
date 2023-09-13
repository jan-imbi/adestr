set.seed(123)
dat <- data.frame(
  endpoint = rnorm(100),
  group = factor(rep(c("ctl", "trt", "ctl", "trt"), c(25,25,25,25))),
  stage = rep(c(1L, 2L), c(50, 50))
)
test_that("Analysis function doesn't throw an error.",
          {
            expect_error(analyze(
              data = dat,
              statistics = list(FirstStageSampleMean(), SampleMean(), NaiveCI()),
              data_distribution = Normal(),
              sigma = 1,
              design = get_example_design()
            ),
            NA)
          })




