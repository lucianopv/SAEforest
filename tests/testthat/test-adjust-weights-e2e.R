# End to end: the weights route through SAEforest_model's nonparametric MSE
# bootstrap, in the combinations the groundtruth var.adjust study runs
# (adjust_mse plugin / refit, re_method insample / crossfit).

weights_e2e_run <- function(...) {
  d <- tiny_saef_data()
  SAEforest_model(Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp,
    pop_data = d$pop, meanOnly = FALSE, MSE = "nonparametric", B = 2,
    num.trees = 50, threshold = median(d$Y), var.adjust = TRUE,
    adjust_method = "weights", B_adj = 5, seed = 1, ...)
}

expect_finite_mse_every_domain <- function(res) {
  d <- tiny_saef_data()
  mse <- res$MSE_Estimates
  expect_setequal(as.character(mse[[d$dName]]), unique(as.character(d$pop[[d$dName]])))
  expect_true(all(is.finite(as.matrix(mse[, -1]))))
  # the weights route really ran: the point fit kept its in-bag counts
  expect_false(is.null(res$MERFmodel$Forest$inbag.counts))
}

test_that("weights route, plug-in K in the MSE bootstrap: finite MSEs for every domain", {
  skip_on_cran(); skip_if_not_installed("emdi")
  expect_finite_mse_every_domain(weights_e2e_run(adjust_mse = "plugin"))
})

test_that("weights route, K refitted in every MSE replicate: finite MSEs for every domain", {
  skip_on_cran(); skip_if_not_installed("emdi")
  expect_finite_mse_every_domain(weights_e2e_run(adjust_mse = "refit"))
})

test_that("weights route with cross-fitted random effects: finite MSEs for every domain", {
  skip_on_cran(); skip_if_not_installed("emdi")
  expect_finite_mse_every_domain(
    weights_e2e_run(adjust_mse = "plugin", re_method = "crossfit", re_folds = 2))
})
