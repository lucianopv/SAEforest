test_that("adj_tol reaches the fitting loop and changes the inner-forest count", {
  d <- tiny_saef_data()
  set.seed(5)
  m_legacy <- SAEforest_model(
    Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp, pop_data = d$pop,
    meanOnly = FALSE, MSE = "none", var.adjust = TRUE, B_adj = 40, adj_tol = 0,
    num.trees = 50, mtry = 3
  )
  set.seed(5)
  m_adapt <- SAEforest_model(
    Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp, pop_data = d$pop,
    meanOnly = FALSE, MSE = "none", var.adjust = TRUE, B_adj = 40, adj_tol = 0.05,
    num.trees = 50, mtry = 3
  )
  # both return valid Indicators; adaptive differs from legacy (fewer inner forests)
  expect_s3_class(m_legacy$Indicators, "data.frame")
  expect_false(isTRUE(all.equal(m_adapt$Indicators, m_legacy$Indicators)))
})

test_that("negative adj_tol is rejected with a clear error", {
  d <- tiny_saef_data()
  expect_error(
    SAEforest_model(
      Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp, pop_data = d$pop,
      meanOnly = FALSE, MSE = "none", var.adjust = TRUE, B_adj = 10, adj_tol = -1,
      num.trees = 50, mtry = 3
    ),
    "adj_tol"
  )
})
