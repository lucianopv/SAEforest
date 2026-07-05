skip_if_serial_only <- function() {
  testthat::skip_on_cran()
  testthat::skip_on_os("windows")
}

base_args <- function(d) list(
  Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp, pop_data = d$pop,
  meanOnly = FALSE, MSE = "nonparametric", B = 2, num.trees = 50, mtry = 3
)

test_that("nonLin MSE is reproducible under a fixed (seed, cores)", {
  skip_if_serial_only()
  d <- tiny_saef_data()
  set.seed(7); r1 <- do.call(SAEforest_model, c(base_args(d), list(cores = 2, seed = 99)))
  set.seed(7); r2 <- do.call(SAEforest_model, c(base_args(d), list(cores = 2, seed = 99)))
  expect_equal(r1$MSE_Estimates, r2$MSE_Estimates)
})

test_that("point estimates are invariant to cores (only MSE differs)", {
  skip_if_serial_only()
  d <- tiny_saef_data()
  set.seed(7); serial <- do.call(SAEforest_model, c(base_args(d), list(cores = 1, seed = 99)))
  set.seed(7); par    <- do.call(SAEforest_model, c(base_args(d), list(cores = 2, seed = 99)))
  expect_equal(par$Indicators, serial$Indicators, tolerance = 1e-8)
})

test_that("cores = 1 equals the plain serial default", {
  d <- tiny_saef_data()
  set.seed(7); a <- do.call(SAEforest_model, c(base_args(d), list(cores = 1, seed = 99)))
  set.seed(7); b <- do.call(SAEforest_model, c(base_args(d), list(seed = 99)))
  expect_equal(a$MSE_Estimates, b$MSE_Estimates, tolerance = 1e-8)
})

test_that("num.threads is overridden with a warning under cores > 1", {
  skip_if_serial_only()
  d <- tiny_saef_data()
  expect_warning(
    do.call(SAEforest_model, c(base_args(d), list(cores = 2, seed = 1, num.threads = 4))),
    "num.threads"
  )
})

test_that("invalid cores is rejected", {
  d <- tiny_saef_data()
  expect_error(do.call(SAEforest_model, c(base_args(d), list(cores = 0))),   "cores")
  expect_error(do.call(SAEforest_model, c(base_args(d), list(cores = 2.5))), "cores")
})
