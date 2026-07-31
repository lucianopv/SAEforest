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

test_that("var.adjust=TRUE rejects B_adj < 1 with a clear error (no NaN crash)", {
  d <- tiny_saef_data()
  expect_error(
    MERFranger(Y = d$Y, X = d$X, random = "(1|district)", data = d$smp,
               var.adjust = TRUE, B_adj = 0, adj_tol = 0.05,
               num.trees = 50, seed = 1),
    "B_adj"
  )
})

# --- Missing OOB predictions in the inner adjustment forests -----------------
# ranger returns NA for an observation that happens to be in-bag in EVERY tree.
# At the fork's default num.trees the per-observation probability is ~1e-10, but
# adjust_ErrorSD_ aggregates an n x B_adj matrix (1.3M cells in the GridSAE
# consumption fit), and neither rowMeans() nor mean() passed na.rm -- so ONE
# missing cell made K NA, and MERFranger's `if (K > naive_unadj^2)` then died
# with "missing value where TRUE/FALSE needed". Observed for real: GridSAE
# src/155, MERF depth5 fold 5, ~19 min into the fit.
#
# num.trees = 1 forces the same condition deterministically: only ~63% of
# observations are in-bag in a single tree, so every inner forest leaves some
# observation without an OOB prediction. The main forest keeps 50 trees, since
# the production NA came from the INNER forests (an NA in rf$predictions would
# instead have cascaded into ranger as a missing-y error, a different failure).

test_that("legacy path keeps K finite when inner forests miss OOB predictions", {
  d <- tiny_saef_data(); rf <- fit_rf(d)
  expect_false(anyNA(rf$predictions))  # the main forest is clean; the inner ones will not be
  set.seed(11)
  res <- adjust_ErrorSD_(Y = d$Y, X = d$X, smp_data = d$smp, rf = rf,
                         B = 4, adj_tol = 0, num.trees = 1)
  expect_true(is.finite(res$K))
  expect_gte(res$K, 0)
})

test_that("adaptive path keeps K finite when inner forests miss OOB predictions", {
  d <- tiny_saef_data(); rf <- fit_rf(d)
  set.seed(12)
  res <- adjust_ErrorSD_(Y = d$Y, X = d$X, smp_data = d$smp, rf = rf,
                         B = 25, adj_tol = 0.05, num.trees = 1)
  expect_true(is.finite(res$K))
  expect_gte(res$K, 0)
})

test_that("a residual pool with no finite entries is reported, not silently averaged", {
  # If EVERY residual is missing, dropping NAs would leave nothing to average and
  # K would come back NaN -- which would resurface as the same if(NA) failure one
  # level up. That case is a real pathology and must fail loudly instead.
  d <- tiny_saef_data(); rf <- fit_rf(d)
  rf_broken <- rf
  rf_broken$predictions <- rep(NA_real_, length(rf$predictions))
  expect_error(
    adjust_ErrorSD_(Y = d$Y, X = d$X, smp_data = d$smp, rf = rf_broken,
                    B = 2, adj_tol = 0, num.trees = 50),
    "no finite OOB residuals"
  )
})
