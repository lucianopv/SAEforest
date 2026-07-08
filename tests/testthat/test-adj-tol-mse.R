test_that("adj_tol_mse defaulting to NULL reproduces the adj_tol MSE result", {
  skip_on_cran(); skip_if_not_installed("emdi")
  d <- tiny_saef_data()
  base <- SAEforest_model(Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp,
    pop_data = d$pop, meanOnly = FALSE, MSE = "nonparametric", B = 2, num.trees = 50,
    threshold = median(d$Y), var.adjust = TRUE, B_adj = 5, adj_tol = 0.05, seed = 1)
  same <- SAEforest_model(Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp,
    pop_data = d$pop, meanOnly = FALSE, MSE = "nonparametric", B = 2, num.trees = 50,
    threshold = median(d$Y), var.adjust = TRUE, B_adj = 5, adj_tol = 0.05,
    adj_tol_mse = NULL, seed = 1)
  expect_equal(base$MSE_Estimates, same$MSE_Estimates)
})

test_that("adj_tol_mse reaches the MSE refits, distinct from adj_tol", {
  skip_on_cran(); skip_if_not_installed("emdi")
  seen <- new.env()
  # point_nonLin is called twice in this run: once for the point estimate (with
  # adj_tol) and once per MSE bootstrap replicate refit (with adj_tol_mse, after
  # resolution). A mock that stops on first capture would only ever see the point
  # estimate's call and never reach the refit. Instead, capture on every call and
  # delegate to the real implementation so the run completes normally; the refit
  # call(s) happen strictly after the point-estimate call, so the last captured
  # value reflects the refit's adj_tol.
  orig_point_nonLin <- point_nonLin
  testthat::local_mocked_bindings(
    point_nonLin = function(..., adj_tol) {
      seen$adj_tol <- adj_tol
      orig_point_nonLin(..., adj_tol = adj_tol)
    }
  )
  d <- tiny_saef_data()
  SAEforest_model(Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp,
    pop_data = d$pop, meanOnly = FALSE, MSE = "nonparametric", B = 2, num.trees = 50,
    threshold = median(d$Y), var.adjust = TRUE, B_adj = 5, adj_tol = 0.05,
    adj_tol_mse = 0.20, seed = 1)
  expect_equal(seen$adj_tol, 0.20)
})
