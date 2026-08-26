# K_fixed: plug-in variance-correction K (spec 2026-08-26, design B).

test_that("MERFranger honours K_fixed and never calls the inner bootstrap", {
  skip_on_cran(); skip_if_not_installed("emdi")
  testthat::local_mocked_bindings(
    adjust_ErrorSD_ = function(...) stop("adjust_ErrorSD_ must not run under K_fixed")
  )
  d <- tiny_saef_data()
  set.seed(7)
  m <- MERFranger(Y = d$Y, X = d$X, random = "(1|district)", data = d$smp,
                  num.trees = 25, var.adjust = TRUE, K_fixed = c(0.4, 0.3))
  expect_gte(m$IterationsUsed, 1L)
  expected <- vapply(seq_len(m$IterationsUsed),
                     function(i) c(0.4, 0.3)[[min(i, 2L)]], numeric(1))
  expect_equal(unlist(m$K), expected)
})

test_that("K_fixed input guards", {
  skip_on_cran(); skip_if_not_installed("emdi")
  d <- tiny_saef_data()
  expect_error(
    MERFranger(Y = d$Y, X = d$X, random = "(1|district)", data = d$smp,
               num.trees = 25, K_fixed = 0.4),
    "requires var.adjust")
  expect_error(
    MERFranger(Y = d$Y, X = d$X, random = "(1|district)", data = d$smp,
               num.trees = 25, var.adjust = TRUE, K_fixed = NA_real_),
    "finite numeric")
  expect_error(
    MERFranger(Y = d$Y, X = d$X, random = "(1|district)", data = d$smp,
               num.trees = 25, var.adjust = TRUE, K_fixed = numeric(0)),
    "finite numeric")
})

plugin_run <- function(am, d = tiny_saef_data(), ...) {
  SAEforest_model(Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp,
    pop_data = d$pop, meanOnly = FALSE, MSE = "nonparametric", B = 2,
    num.trees = 50, threshold = median(d$Y), var.adjust = TRUE, B_adj = 5,
    adj_tol = 0.05, seed = 1, adjust_mse = am, ...)
}

test_that("adjust_mse='plugin': same point estimates, finite MSEs, seeded-reproducible", {
  skip_on_cran(); skip_if_not_installed("emdi")
  refit  <- plugin_run("refit")
  plugin <- plugin_run("plugin")
  # the point fit is shared and runs before any plugin divergence
  expect_identical(refit$Indicators, plugin$Indicators)
  expect_true(all(is.finite(as.matrix(plugin$MSE_Estimates[, -1]))))
  # plugin consumes no RNG in the adjustment, so streams (and MSEs) differ
  expect_false(identical(refit$MSE_Estimates, plugin$MSE_Estimates))
  # but the plugin estimator itself is exactly reproducible under the seed
  expect_identical(plugin$MSE_Estimates, plugin_run("plugin")$MSE_Estimates)
})

test_that("the point fit's K values reach the replicate fits", {
  skip_on_cran(); skip_if_not_installed("emdi")
  seen <- new.env()
  orig_point_nonLin <- point_nonLin
  testthat::local_mocked_bindings(
    point_nonLin = function(..., K_fixed = NULL) {
      # point-estimate call sees NULL; replicate calls run strictly after it,
      # so the last captured value is the replicates' K_fixed.
      seen$last_K <- K_fixed
      orig_point_nonLin(..., K_fixed = K_fixed)
    }
  )
  res <- plugin_run("plugin")
  expect_true(is.numeric(seen$last_K) && length(seen$last_K) >= 1)
  expect_true(all(is.finite(seen$last_K)))
  # The replicates' K_fixed is exactly the point fit's own per-iteration K
  # (the mock delegates, so the run returns normally; MERFmodel is a c() of
  # the MERFranger object with extras, so $K lives there unchanged).
  expect_equal(seen$last_K, unlist(res$MERFmodel$K))
})

test_that("sanitize_K_fixed substitutes non-finite point-fit K with 0, loudly", {
  expect_warning(
    out <- sanitize_K_fixed(c(0.4, NA_real_, 0.3)),
    "Non-finite point-fit K"
  )
  expect_equal(out, c(0.4, 0, 0.3))
  expect_silent(out2 <- sanitize_K_fixed(c(0.4, 0.3)))
  expect_identical(out2, c(0.4, 0.3))
})

test_that("adjust_mse guards", {
  skip_on_cran(); skip_if_not_installed("emdi")
  d <- tiny_saef_data()
  expect_error(
    SAEforest_model(Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp,
      pop_data = d$pop, meanOnly = FALSE, MSE = "nonparametric", B = 2,
      num.trees = 25, threshold = median(d$Y), var.adjust = FALSE, seed = 1,
      adjust_mse = "plugin"),
    "requires var.adjust")
  expect_error(
    SAEforest_model(Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp,
      pop_data = d$pop, meanOnly = FALSE, MSE = "none", num.trees = 25,
      threshold = median(d$Y), var.adjust = TRUE, seed = 1,
      adjust_mse = "plugin"),
    "MSE")
})
