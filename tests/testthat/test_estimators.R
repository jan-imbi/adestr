# Test everything that hasn't been tested somehwere else here.
test_that("MidpointEstimators don't error",
{expect_error(
  {
    ests <- get_stagewise_estimators(MidpointMLEOrderingCI(), Normal(FALSE), FALSE, designad, 1, FALSE)
    ests[[1]](smean1= .1, sigma=1, two_armed=FALSE, n1=designad@n1, design=designad)
    ests[[2]](smean1= .1, smean2=.2, sigma=1, two_armed=FALSE, n1=designad@n1, n2=30, design=designad)
  }
  , NA
  )
  expect_error(
    {
      ests <- get_stagewise_estimators(MidpointLikelihoodRatioOrderingCI(), Normal(FALSE), FALSE, designad, 1, FALSE)
      ests[[1]](smean1= .1, sigma=1, two_armed=FALSE, n1=designad@n1, design=designad)
      ests[[2]](smean1= .1, smean2=.2, sigma=1, two_armed=FALSE, n1=designad@n1, n2=30, design=designad)
    }
    , NA
  )
  expect_error(
    {
      ests <- get_stagewise_estimators(MidpointScoreTestOrderingCI(), Normal(FALSE), FALSE, designad, 1, FALSE)
      ests[[1]](smean1= .1, sigma=1, two_armed=FALSE, n1=designad@n1, design=designad)
      ests[[2]](smean1= .1, smean2=.2, sigma=1, two_armed=FALSE, n1=designad@n1, n2=30, design=designad)
    }
    , NA
  )
  expect_error(
    {
      ests <- get_stagewise_estimators(MidpointStagewiseCombinationFunctionOrderingCI(), Normal(FALSE), FALSE, designad, 1, FALSE)
      ests[[1]](smean1= .1, sigma=1, two_armed=FALSE, n1=designad@n1, design=designad)
      ests[[2]](smean1= .1, smean2=.2, sigma=1, two_armed=FALSE, n1=designad@n1, n2=30, design=designad)
    }
    , NA
  )
})


