test_that("point + MSE output has one row per domain and valid indicator ranges", {
  d <- tiny_saef_data()
  res <- SAEforest_model(
    Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp, pop_data = d$pop,
    meanOnly = FALSE, MSE = "nonparametric", var.adjust = TRUE,
    B = 2, B_adj = 2, num.trees = 50, mtry = 3, seed = 6
  )
  n_dom <- length(unique(d$pop$district))
  expect_equal(nrow(res$Indicators), n_dom)
  expect_equal(nrow(res$MSE_Estimates), n_dom)
  expect_true(all(res$Indicators$Hcr >= 0 & res$Indicators$Hcr <= 1))
  expect_true(all(res$Indicators$Gini >= 0 & res$Indicators$Gini <= 1))
  mse_num <- res$MSE_Estimates[, sapply(res$MSE_Estimates, is.numeric)]
  expect_false(anyNA(as.matrix(mse_num)))
  expect_true(all(as.matrix(mse_num) >= 0))
})
