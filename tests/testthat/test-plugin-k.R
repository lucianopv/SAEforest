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
