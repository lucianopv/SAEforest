# Tests for the unsupported-combo warnings added in R/SAEforest_model.R:
# var.adjust/transformation are only implemented on the non-linear smearing path
# and are silently ignored elsewhere unless explicitly warned about.

test_that("var.adjust = TRUE with meanOnly = TRUE warns that it is ignored", {
  d <- tiny_saef_data()
  expect_warning(
    SAEforest_model(
      Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp, pop_data = d$pop,
      meanOnly = TRUE, MSE = "none", var.adjust = TRUE,
      num.trees = 50, mtry = 3, seed = 1
    ),
    "ignored"
  )
})

test_that("transformation = 'log' with meanOnly = TRUE warns that it is ignored", {
  d <- tiny_saef_data()
  expect_warning(
    SAEforest_model(
      Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp, pop_data = d$pop,
      meanOnly = TRUE, MSE = "none", transformation = "log",
      num.trees = 50, mtry = 3, seed = 1
    ),
    "ignored"
  )
})

test_that("transformation = 'log' with smearing = FALSE warns that it is ignored", {
  d <- tiny_saef_data()
  expect_warning(
    SAEforest_model(
      Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp, pop_data = d$pop,
      meanOnly = FALSE, MSE = "none", transformation = "log", smearing = FALSE,
      num.trees = 50, mtry = 3, seed = 1
    ),
    "ignored"
  )
})

test_that("a normal non-linear var.adjust = TRUE call raises no unsupported-combo warning", {
  d <- tiny_saef_data()
  expect_no_warning_matching <- function(expr, pattern) {
    warnings_seen <- character(0)
    withCallingHandlers(
      expr,
      warning = function(w) {
        warnings_seen[[length(warnings_seen) + 1]] <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      }
    )
    expect_false(any(grepl(pattern, warnings_seen)))
  }

  expect_no_warning_matching(
    SAEforest_model(
      Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp, pop_data = d$pop,
      meanOnly = FALSE, MSE = "none", var.adjust = TRUE,
      num.trees = 50, mtry = 3, seed = 1
    ),
    "ignored"
  )
})
