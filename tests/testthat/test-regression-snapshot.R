test_that("adjusted point estimates match the stored golden reference", {
  d <- tiny_saef_data()
  ref <- readRDS(testthat::test_path("fixtures", "adj_nonLin_none.rds"))
  set.seed(2026)
  got <- SAEforest_model(
    Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp, pop_data = d$pop,
    meanOnly = FALSE, MSE = "none", var.adjust = TRUE, B_adj = 2,
    num.trees = 50, mtry = 3
  )$Indicators
  expect_equal(got, ref, tolerance = 1e-8)
})

test_that("legacy path (adj_tol = 0) matches the frozen pre-adaptive reference", {
  d <- tiny_saef_data()
  ref <- readRDS(testthat::test_path("fixtures", "adj_nonLin_none_legacy.rds"))
  set.seed(2026)
  got <- SAEforest_model(
    Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp, pop_data = d$pop,
    meanOnly = FALSE, MSE = "none", var.adjust = TRUE, B_adj = 2, adj_tol = 0,
    num.trees = 50, mtry = 3
  )$Indicators
  expect_equal(got, ref, tolerance = 1e-8)
})
