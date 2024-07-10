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

test_that("Analysis function throws an error when it should.",
          {
            expect_warning(
              analyze(
                data = dat,
                statistics = list(FirstStageSampleMean(), SampleMean(), NaiveCI()),
                data_distribution = Student(TRUE),
                sigma = 1,
                design = get_example_design(TRUE)
              ))
            expect_warning(expect_warning(expect_warning(
              analyze(
                data = dat,
                statistics = list(FirstStageSampleMean(), SampleMean(), NaiveCI()),
                data_distribution = Normal(TRUE),
                design = get_example_design(TRUE)
              ))))
            expect_error(analyze(
              data = dat,
              statistics = list(FirstStageSampleMean(), SampleMean(), NaiveCI()),
              data_distribution = Student(FALSE),
              design = get_example_design(TRUE)
            ))

          })

test_that("Single-stage analyze works",
{
  expect_warning(
    expect_error(analyze(
      data = dat[dat$stage==1,],
      statistics = list(FirstStageSampleMean(), SampleMean(), NaiveCI()),
      data_distribution = Normal(TRUE),
      sigma = 1,
      design = get_example_design(TRUE)
    ),
    NA)
    )
})


