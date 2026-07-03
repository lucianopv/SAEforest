test_that("var.adjust=FALSE reproduces the stored upstream reference", {
  d <- tiny_saef_data()
  ref <- readRDS(testthat::test_path("fixtures", "upstream_nonLin_none.rds"))
  set.seed(2026)
  got <- SAEforest_model(
    Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp, pop_data = d$pop,
    meanOnly = FALSE, MSE = "none", num.trees = 50, mtry = 3
  )$Indicators
  expect_equal(got, ref, tolerance = 1e-8)
})
