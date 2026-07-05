# End-to-end coverage for the select.indicator speed option: it must flow through
# point estimation AND the MSE bootstrap, and an invalid name must fail early with
# a clear message (not a cryptic subscript error deep inside the bootstrap).

test_that("select.indicator returns only the requested indicator end-to-end (point + MSE)", {
  d <- tiny_saef_data()
  set.seed(6)
  res <- SAEforest_model(
    Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp, pop_data = d$pop,
    meanOnly = FALSE, MSE = "nonparametric", select.indicator = "Hcr",
    B = 2, B_adj = 2, num.trees = 50, mtry = 3
  )

  # Indicators: the domain id column plus ONLY the selected indicator.
  expect_true("Hcr" %in% names(res$Indicators))
  expect_false(any(c("Gini", "Mean", "Qsr", "Median") %in% names(res$Indicators)))
  expect_equal(nrow(res$Indicators), length(unique(d$pop$district)))
  expect_true(all(res$Indicators$Hcr >= 0 & res$Indicators$Hcr <= 1))

  # MSE estimates are restricted to the same indicator.
  expect_true("Hcr" %in% names(res$MSE_Estimates))
  expect_false(any(c("Gini", "Mean") %in% names(res$MSE_Estimates)))
})

test_that("an invalid select.indicator errors early with a clear message", {
  d <- tiny_saef_data()
  # MSE = "none" so no bootstrap runs: the error must surface during point
  # estimation, not after an expensive bootstrap.
  expect_error(
    SAEforest_model(
      Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp, pop_data = d$pop,
      meanOnly = FALSE, MSE = "none", select.indicator = "Povrty",
      num.trees = 50, mtry = 3
    ),
    "select.indicator"
  )
})
