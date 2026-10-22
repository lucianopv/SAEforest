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

test_that("NA out-of-bag predictions keep every row on its own area's random effect", {
  # At 6 trees about 6% of rows are in-bag in every tree, so their OOB
  # prediction, r_ij and r_ij_new are NA. lmer drops those rows, so
  # predict(lmefit) was shorter than the data and Target - AllEffects recycled
  # (rows after the first NA took another row's effect), and so did the pool.
  d <- make_toy()
  seen <- new.env(); seen$Y <- list(); seen$n_na <- integer(0)
  orig <- adjust_ErrorSD_weights_
  testthat::local_mocked_bindings(adjust_ErrorSD_weights_ = function(Y, X, rf, ...) {
    seen$Y[[length(seen$Y) + 1L]] <- Y          # AdjustedTarget of this iteration
    seen$n_na <- c(seen$n_na, sum(!is.finite(rf$predictions)))
    orig(Y, X, rf, ...)
  })
  out <- collect_warnings(
    MERFranger(Y = d$Y, X = d$X, random = "(1|area)", data = d$data,
               var.adjust = TRUE, adjust_method = "weights", B_adj = 10,
               num.trees = 6, MaxIterations = 4, num.threads = 1, seed = 3))
  m <- out$value
  expect_true(all(seen$n_na > 0))                # the case is present every iteration
  expect_length(out$warnings, 0)                 # no "longer object length" recycling
  # Each iteration after the first is given Target - AllEffects: one value per
  # row, the same within an area (that area's random effect).
  for (Y in seen$Y[-1]) {
    expect_length(Y, length(d$Y))
    spread <- tapply(d$Y - Y, d$data$area, function(e) diff(range(e)))
    expect_lt(max(spread), 1e-10)
  }
  # The final pool: each fitted row's r_ij_new less its own area's random effect.
  fit_rows <- as.integer(rownames(stats::model.frame(m$EffectModel)))
  re <- m$RandomEffects$area[as.character(d$data$area[fit_rows]), 1]
  pool <- lme4::getME(m$EffectModel, "y") - re
  expect_equal(unname(m$OOBresiduals), unname(pool - mean(pool)), tolerance = 1e-10)
})
