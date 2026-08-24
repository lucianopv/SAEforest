# Byte-parity guards for the numerically inert MSE-bootstrap edits:
# the `groups` hoist (Task 3) and the pseudo-population memory edits (Task 4).
# Every change these tests cover must leave the MSE estimates untouched.

mse_nonlin_run <- function(seed = 2026) {
  d <- tiny_saef_data()
  SAEforest_model(
    Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp, pop_data = d$pop,
    meanOnly = FALSE, MSE = "nonparametric", B = 3,
    threshold = median(d$Y), var.adjust = TRUE, B_adj = 5, adj_tol = 0.05,
    num.trees = 50, mtry = 3, seed = seed
  )$MSE_Estimates
}

mse_mean_run <- function(seed = 2026) {
  d <- tiny_saef_data()
  SAEforest_model(
    Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp, pop_data = d$pop,
    meanOnly = TRUE, MSE = "nonparametric", B = 3,
    num.trees = 50, mtry = 3, seed = seed
  )$MSE_Estimates
}

# The aggregate_to branch of my_agg (MSE_SAEforest_nonLin.R) -- dName finer
# than aggregate_to, e.g. dName="subunit"/aggregate_to="district". This is the
# only branch every GridSAE fit exercises; mse_nonlin_run()/mse_mean_run()
# above only cover aggregate_to = NULL. Added per final review 2026-08-24
# (finding M2): the Map(cbind, ...) removal rewrote both branches of my_agg,
# but nothing in the suite previously combined MSE = "nonparametric" with
# aggregate_to != NULL.
mse_nonlin_agg_run <- function(seed = 2026) {
  d <- tiny_saef_data_with_subunit()
  SAEforest_model(
    Y = d$Y, X = d$X, dName = "subunit", smp_data = d$smp, pop_data = d$pop,
    meanOnly = FALSE, MSE = "nonparametric", B = 3,
    threshold = median(d$Y), var.adjust = TRUE, B_adj = 5, adj_tol = 0.05,
    num.trees = 50, mtry = 3, seed = seed,
    aggregate_to = "district", select.indicator = c("Mean", "Hcr")
  )$MSE_Estimates
}

test_that("nonlinear MSE bootstrap is reproducible for a fixed seed", {
  skip_on_cran(); skip_if_not_installed("emdi")
  expect_equal(mse_nonlin_run(seed = 5), mse_nonlin_run(seed = 5))
})

test_that("mean MSE bootstrap is reproducible for a fixed seed", {
  skip_on_cran(); skip_if_not_installed("emdi")
  expect_equal(mse_mean_run(seed = 5), mse_mean_run(seed = 5))
})

test_that("nonlinear MSE bootstrap matches the frozen pre-hoist golden", {
  skip_on_cran(); skip_if_not_installed("emdi")
  ref <- readRDS(test_path("fixtures", "mse_memory_parity_ref.rds"))
  expect_equal(mse_nonlin_run(seed = 5), ref$nonlin, tolerance = 0)
})

test_that("mean MSE bootstrap matches the frozen pre-hoist golden", {
  skip_on_cran(); skip_if_not_installed("emdi")
  ref <- readRDS(test_path("fixtures", "mse_memory_parity_ref.rds"))
  expect_equal(mse_mean_run(seed = 5), ref$mean, tolerance = 0)
})

test_that("nonlinear MSE bootstrap with aggregate_to matches the frozen golden", {
  skip_on_cran(); skip_if_not_installed("emdi")
  ref <- readRDS(test_path("fixtures", "mse_memory_parity_ref.rds"))
  expect_equal(mse_nonlin_agg_run(seed = 5), ref$nonlin_agg, tolerance = 0)
})
