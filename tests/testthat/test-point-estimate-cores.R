test_that("point_nonLin produces identical Indicators for cores = 1 and cores = 2", {
  skip_on_os("windows")
  d <- tiny_saef_data()

  set.seed(21)
  serial <- point_nonLin(
    Y = d$Y, X = d$X, dName = d$dName, threshold = 15000,
    smp_data = d$smp, pop_data = d$pop,
    initialRandomEffects = 0, ErrorTolerance = 0.0001, MaxIterations = 25,
    custom_indicator = NULL, aggregate_to = NULL,
    transformation = "log", num.trees = 50, cores = 1
  )

  set.seed(21)
  parallel_result <- point_nonLin(
    Y = d$Y, X = d$X, dName = d$dName, threshold = 15000,
    smp_data = d$smp, pop_data = d$pop,
    initialRandomEffects = 0, ErrorTolerance = 0.0001, MaxIterations = 25,
    custom_indicator = NULL, aggregate_to = NULL,
    transformation = "log", num.trees = 50, cores = 2
  )

  expect_equal(serial[[1]], parallel_result[[1]])
})
