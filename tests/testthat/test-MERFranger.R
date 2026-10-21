# Test if MERFranger and most important errors work as expected

test_that("MERFranger and prediction generics work?", {
  data("eusilcA_pop")
  data("eusilcA_smp")

  income <- eusilcA_smp$eqIncome
  X_covar <- eusilcA_smp[, -c(1, 16, 17, 18)]

  model1 <- MERFranger(
    Y = income, X = X_covar, random = "(1|district)",
    data = eusilcA_smp, num.trees = 50
  )

  # get individual predictions:
  ind_pred <- predict(model1, eusilcA_pop)
  ind_pred_no <- predict(model1$Forest, eusilcA_pop)$predictions + predict(model1$EffectModel,
    eusilcA_pop,
    allow.new.levels = TRUE
  )

  expect_s3_class(model1, "MERFranger")
  expect_equal(ind_pred, ind_pred_no)
})


test_that("MERFranger produces informative error messages?", {
  data("eusilcA_pop")
  data("eusilcA_smp")

  income <- 0
  X_covar <- eusilcA_smp[, -c(1, 16, 17, 18)]

  expect_error(MERFranger(
    Y = income, X = X_covar, random = "(1|district)",
    data = eusilcA_smp, num.trees = 50
  ))

  income <- eusilcA_smp$eqIncome[1:100]
  X_covar <- eusilcA_smp[, -c(1, 16, 17, 18)]

  expect_error(MERFranger(
    Y = income, X = X_covar, random = "(1|district)",
    data = eusilcA_smp, num.trees = 50
  ))
})

# var.adjust EM loop: degenerate cases (2026-10-22 review of the weights route).

# Collects warnings instead of letting them escape the test.
collect_warnings <- function(expr) {
  w <- character(0)
  value <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  })
  list(value = value, warnings = w)
}

test_that("a clamped EM iteration ends the loop instead of crashing it", {
  # K far above the naive residual variance (about 2) clamps the adjusted SD to
  # 0, so lmer sees an all-zero response and returns logLik = +Inf. The next
  # convergence check divided by that Inf and `while (NA)` threw.
  d <- make_toy()
  out <- collect_warnings(
    MERFranger(Y = d$Y, X = d$X, random = "(1|area)", data = d$data,
               var.adjust = TRUE, K_fixed = 50, num.trees = 50,
               num.threads = 1, seed = 3))
  m <- out$value
  expect_s3_class(m, "MERFranger")
  expect_equal(m$ErrorSD, 0)
  expect_equal(m$RanEffSD, 0)
  expect_true(is.integer(m$clamped) && m$clamped >= 1L)
  # one clamp warning per clamped iteration, and nothing else
  expect_true(all(grepl("^Variance bias correction K exceeds naive residual variance; clamping adjusted SD to 0",
                        out$warnings)))
  expect_equal(length(out$warnings), m$clamped)
})

test_that("clamped is 0 on a fit that never clamps, and absent without var.adjust", {
  d <- make_toy()
  common <- list(Y = d$Y, X = d$X, random = "(1|area)", data = d$data,
                 num.trees = 50, num.threads = 1, seed = 3)
  m_adj <- do.call(MERFranger, c(common, list(var.adjust = TRUE, K_fixed = 0.3)))
  expect_identical(m_adj$clamped, 0L)
  expect_null(do.call(MERFranger, common)$clamped)
})

test_that("non-degenerate var.adjust fits are unchanged by the degenerate-case guards", {
  ref <- readRDS(testthat::test_path("fixtures", "varadjust_fits_ref.rds"))
  expect_identical(varadjust_guard_fits(), ref)
})
