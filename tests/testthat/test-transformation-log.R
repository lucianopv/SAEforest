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
  expect_error(
    SAEforest_model(
      Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp, pop_data = d$pop,
      meanOnly = FALSE, MSE = "nonparametric", transformation = "log",
      B = 2, B_adj = 2, num.trees = 50, mtry = 3, seed = 4
    ),
    NA
  )
})
