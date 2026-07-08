test_that("n_smear_residuals = NULL (default) reproduces the exact unmodified smearing result", {
  d <- tiny_saef_data()
  set.seed(31)
  default_result <- point_nonLin(
    Y = d$Y, X = d$X, dName = d$dName, threshold = 15000,
    smp_data = d$smp, pop_data = d$pop,
    initialRandomEffects = 0, ErrorTolerance = 0.0001, MaxIterations = 25,
    custom_indicator = NULL, aggregate_to = NULL,
    transformation = "log", num.trees = 50, n_smear_residuals = NULL
  )
  set.seed(31)
  explicit_null_result <- point_nonLin(
    Y = d$Y, X = d$X, dName = d$dName, threshold = 15000,
    smp_data = d$smp, pop_data = d$pop,
    initialRandomEffects = 0, ErrorTolerance = 0.0001, MaxIterations = 25,
    custom_indicator = NULL, aggregate_to = NULL,
    transformation = "log", num.trees = 50
  )
  expect_equal(default_result[[1]], explicit_null_result[[1]])
})

test_that("n_smear_residuals subsamples reproducibly under a fixed seed", {
  d <- tiny_saef_data()

  set.seed(32)
  run1 <- point_nonLin(
    Y = d$Y, X = d$X, dName = d$dName, threshold = 15000,
    smp_data = d$smp, pop_data = d$pop,
    initialRandomEffects = 0, ErrorTolerance = 0.0001, MaxIterations = 25,
    custom_indicator = NULL, aggregate_to = NULL,
    transformation = "log", num.trees = 50, n_smear_residuals = 10
  )

  set.seed(32)
  run2 <- point_nonLin(
    Y = d$Y, X = d$X, dName = d$dName, threshold = 15000,
    smp_data = d$smp, pop_data = d$pop,
    initialRandomEffects = 0, ErrorTolerance = 0.0001, MaxIterations = 25,
    custom_indicator = NULL, aggregate_to = NULL,
    transformation = "log", num.trees = 50, n_smear_residuals = 10
  )

  expect_equal(run1[[1]], run2[[1]])
})

test_that("n_smear_residuals >= the actual residual count is a no-op (uses all residuals, matches default)", {
  d <- tiny_saef_data()
  n_residuals <- nrow(d$smp)

  set.seed(33)
  default_result <- point_nonLin(
    Y = d$Y, X = d$X, dName = d$dName, threshold = 15000,
    smp_data = d$smp, pop_data = d$pop,
    initialRandomEffects = 0, ErrorTolerance = 0.0001, MaxIterations = 25,
    custom_indicator = NULL, aggregate_to = NULL,
    transformation = "log", num.trees = 50, n_smear_residuals = NULL
  )
  set.seed(33)
  capped_result <- point_nonLin(
    Y = d$Y, X = d$X, dName = d$dName, threshold = 15000,
    smp_data = d$smp, pop_data = d$pop,
    initialRandomEffects = 0, ErrorTolerance = 0.0001, MaxIterations = 25,
    custom_indicator = NULL, aggregate_to = NULL,
    transformation = "log", num.trees = 50, n_smear_residuals = n_residuals + 1000
  )
  expect_equal(default_result[[1]], capped_result[[1]])
})
