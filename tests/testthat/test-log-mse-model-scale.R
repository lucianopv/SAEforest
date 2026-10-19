# transformation = "log" fits the forest to log(Y); the MSE bootstrap's error-SD
# correction and its residual pools must take the outcome on that same scale.
# (Before 1.0.1.9002 they took the level-scale Y: the bias term swamped the
# log-scale ErrorSD^2, the adjusted SD was clamped to 0 and the bootstrap drew
# no household errors.)

log_mse_fit <- function(d, Y, transformation, MSE = "nonparametric") {
  suppressMessages(SAEforest_model(
    Y = Y, X = d$X, dName = d$dName, smp_data = d$smp, pop_data = d$pop,
    meanOnly = FALSE, MSE = MSE, transformation = transformation,
    B = 2, B_adj = 20, num.trees = 50, mtry = 3, seed = 4
  ))
}

test_that("log MSE: the adjusted error SD is positive and not clamped", {
  skip_on_cran()
  d <- tiny_saef_data(n_dom = 6, per_dom = 25)
  warns <- character()
  fit <- withCallingHandlers(log_mse_fit(d, d$Y, "log"), warning = function(w) {
    warns <<- c(warns, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_false(any(grepl("clamping", warns)))
  expect_gt(fit$AdjustedSD, 0)
})

test_that("log MSE: the adjusted error SD equals the one from fitting log(Y) untransformed", {
  skip_on_cran()
  d <- tiny_saef_data(n_dom = 6, per_dom = 25)
  # SAEforest_model() wants Y as a column of smp_data
  d$smp$log_eqIncome <- log(d$smp$eqIncome)
  fit_log <- log_mse_fit(d, d$smp$eqIncome, "log")
  fit_none <- log_mse_fit(d, d$smp$log_eqIncome, "none")
  # same point fit on the same numbers, so the same correction
  expect_equal(fit_log$MERFmodel$ErrorSD, fit_none$MERFmodel$ErrorSD)
  expect_equal(fit_log$AdjustedSD, fit_none$AdjustedSD)
})

test_that("log MSE: the bootstrap draws household errors from model-scale residuals", {
  skip_on_cran()
  d <- tiny_saef_data(n_dom = 6, per_dom = 25)
  seen <- NULL
  real_ran_comp <- ran_comp
  testthat::local_mocked_bindings(
    ran_comp = function(mod, smp_data, Y, dName, ADJsd) {
      out <- real_ran_comp(mod = mod, smp_data = smp_data, Y = Y, dName = dName, ADJsd = ADJsd)
      if (is.null(seen)) seen <<- list(Y = Y, forest_res = out$forest_res, ADJsd = ADJsd)
      out
    }
  )
  log_mse_fit(d, d$Y, "log")
  expect_equal(sort(seen$Y), sort(log(d$Y)))
  expect_gt(sd(seen$forest_res), 0)
  expect_equal(sd(seen$forest_res), seen$ADJsd)
})

test_that("log wild MSE: the residuals are taken on the model's scale", {
  skip_on_cran()
  d <- tiny_saef_data(n_dom = 6, per_dom = 25)
  # stop at ran_comp_wild(): only the outcome it receives is checked here
  testthat::local_mocked_bindings(
    ran_comp_wild = function(mod, smp_data, Y, dName, ADJsd) {
      stop(structure(class = c("captured_y", "error", "condition"),
                     list(message = "captured", call = NULL, Y = Y)))
    }
  )
  err <- expect_error(suppressWarnings(log_mse_fit(d, d$Y, "log", MSE = "wild")),
                      class = "captured_y")
  expect_equal(sort(err$Y), sort(log(d$Y)))
})
