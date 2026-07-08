mse_run <- function(mse_tol = 0, seed = 1) {
  d <- tiny_saef_data()
  SAEforest_model(Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp,
    pop_data = d$pop, meanOnly = FALSE, MSE = "nonparametric", B = 4, num.trees = 50,
    threshold = median(d$Y), var.adjust = TRUE, B_adj = 5, adj_tol = 0.05,
    mtry = 3, seed = seed, mse_tol = mse_tol)$MSE_Estimates
}

test_that("mse_tol = 0 reproduces the frozen full-B golden", {
  skip_on_cran(); skip_if_not_installed("emdi")
  golden <- readRDS(test_path("fixtures", "mse_nonlin_default.rds"))
  expect_equal(mse_run(mse_tol = 0), golden, tolerance = 1e-8)
})

test_that("mse_tol > 0 is reproducible for a fixed seed", {
  skip_on_cran(); skip_if_not_installed("emdi")
  expect_equal(mse_run(mse_tol = 0.2, seed = 7), mse_run(mse_tol = 0.2, seed = 7))
})

test_that("B = 2 does not error under mse_tol > 0", {
  skip_on_cran(); skip_if_not_installed("emdi")
  d <- tiny_saef_data()
  expect_no_error(SAEforest_model(Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp,
    pop_data = d$pop, meanOnly = FALSE, MSE = "nonparametric", B = 2, num.trees = 50,
    threshold = median(d$Y), var.adjust = TRUE, B_adj = 5, mtry = 3, seed = 1, mse_tol = 0.2))
})
