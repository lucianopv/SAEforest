test_that("mean-path MSE bootstrap honors the caller's forest config (... forwarded)", {
  d <- tiny_saef_data()
  ref <- readRDS(testthat::test_path("fixtures", "mean_mse_none.rds"))
  set.seed(2026)
  got <- suppressMessages(SAEforest_model(
    Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp, pop_data = d$pop,
    meanOnly = TRUE, MSE = "nonparametric", B = 2, num.trees = 50, mtry = 3
  ))$MSE_Estimates
  expect_equal(got, ref, tolerance = 1e-8)
})
