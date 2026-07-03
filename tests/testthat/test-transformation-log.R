test_that("transformation='log' runs and returns positive, finite indicators", {
  d <- tiny_saef_data()
  out <- SAEforest_model(
    Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp, pop_data = d$pop,
    meanOnly = FALSE, MSE = "none", transformation = "log",
    num.trees = 50, mtry = 3, seed = 4
  )$Indicators
  num <- out[, sapply(out, is.numeric)]
  expect_true(all(is.finite(as.matrix(num))))
  expect_true(all(num[, "Mean"] > 0))
})

test_that("log MSE path (the old double-comma bug) executes without error", {
  d <- tiny_saef_data()
  # On this tiny fixture with B_adj = 2, adjust_ErrorSD()'s bootstrap bias estimate
  # is noisy and legitimately exceeds ErrorSD^2, triggering the expected
  # "clamping adjusted SD to 0" guard (see R/adjust_ErrorSD.R). That guard is
  # benign and orthogonal to what this test checks (absence of an error), so it
  # is suppressed here to keep test output pristine.
  expect_error(
    suppressWarnings(
      SAEforest_model(
        Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp, pop_data = d$pop,
        meanOnly = FALSE, MSE = "nonparametric", transformation = "log",
        B = 2, B_adj = 2, num.trees = 50, mtry = 3, seed = 4
      )
    ),
    NA
  )
})
