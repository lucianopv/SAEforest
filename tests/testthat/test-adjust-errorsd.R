# Fit a small forest on the tiny fixture to feed adjust_ErrorSD_ directly.
fit_rf <- function(d, seed = 7) {
  set.seed(seed)
  ranger::ranger(y = d$Y, x = d$X, num.trees = 50)
}

test_that("adj_tol = 0 fits exactly B forests and returns a list(K, m)", {
  d <- tiny_saef_data(); rf <- fit_rf(d)
  set.seed(1)
  res <- adjust_ErrorSD_(Y = d$Y, X = d$X, smp_data = d$smp, rf = rf,
                         B = 8, adj_tol = 0, num.trees = 50)
  expect_type(res, "list")
  expect_named(res, c("K", "m"))
  expect_equal(res$m, 8)
  expect_true(is.finite(res$K) && res$K >= 0)
})

test_that("adaptive path stops between min(B_MIN,B) and B, K stays close to full-B K", {
  d <- tiny_saef_data(); rf <- fit_rf(d)
  set.seed(2)
  full <- adjust_ErrorSD_(Y = d$Y, X = d$X, smp_data = d$smp, rf = rf,
                          B = 50, adj_tol = 0, num.trees = 50)
  set.seed(2)
  adap <- adjust_ErrorSD_(Y = d$Y, X = d$X, smp_data = d$smp, rf = rf,
                          B = 50, adj_tol = 0.05, num.trees = 50)
  expect_gte(adap$m, min(20, 50))
  expect_lte(adap$m, 50)
  # adaptive K within a loose multiple of adj_tol of the full-cap K
  expect_lt(abs(adap$K - full$K) / full$K, 0.25)
})

test_that("adaptive path is reproducible under a fixed seed", {
  d <- tiny_saef_data(); rf <- fit_rf(d)
  set.seed(3); a <- adjust_ErrorSD_(d$Y, d$X, d$smp, rf, B = 50, adj_tol = 0.05, num.trees = 50)
  set.seed(3); b <- adjust_ErrorSD_(d$Y, d$X, d$smp, rf, B = 50, adj_tol = 0.05, num.trees = 50)
  expect_equal(a, b)
})

test_that("NULL adj_tol behaves as 0 (legacy)", {
  d <- tiny_saef_data(); rf <- fit_rf(d)
  set.seed(4); res <- adjust_ErrorSD_(d$Y, d$X, d$smp, rf, B = 6, adj_tol = NULL, num.trees = 50)
  expect_equal(res$m, 6)
})

test_that("adaptive path survives degenerate B = 1 (sd of length-1 g is NA)", {
  d <- tiny_saef_data(); rf <- fit_rf(d)
  set.seed(5)
  res <- adjust_ErrorSD_(d$Y, d$X, d$smp, rf, B = 1, adj_tol = 0.05, num.trees = 50)
  expect_equal(res$m, 1)
  expect_true(is.finite(res$K) && res$K >= 0)
})

test_that("MERFranger runs with var.adjust=TRUE, B_adj=1, adj_tol>0 (no crash)", {
  d <- tiny_saef_data()
  expect_no_error(
    mod <- MERFranger(Y = d$Y, X = d$X, random = "(1|district)", data = d$smp,
                      var.adjust = TRUE, B_adj = 1, adj_tol = 0.05,
                      num.trees = 50, seed = 1)
  )
  expect_true(is.finite(mod$ErrorSD))
})
