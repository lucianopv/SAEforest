test_that("var.adjust=TRUE uses the in-loop corrected ErrorSD (no MSE double-correction)", {
  d <- tiny_saef_data()
  set.seed(6)
  res <- SAEforest_model(
    Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp, pop_data = d$pop,
    meanOnly = FALSE, MSE = "nonparametric", var.adjust = TRUE,
    B = 2, B_adj = 2, num.trees = 50, mtry = 3
  )
  expect_equal(res$AdjustedSD, res$MERFmodel$ErrorSD, tolerance = 1e-8)
})

test_that("var.adjust=FALSE still bias-corrects the ErrorSD at the MSE stage", {
  d <- tiny_saef_data()
  set.seed(6)
  res <- SAEforest_model(
    Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp, pop_data = d$pop,
    meanOnly = FALSE, MSE = "nonparametric", var.adjust = FALSE,
    B = 2, B_adj = 2, num.trees = 50, mtry = 3
  )
  expect_lte(res$AdjustedSD, res$MERFmodel$ErrorSD)
})
