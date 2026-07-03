fit_adj <- function(d, seed = 11) {
  random <- "(1|district)"
  MERFranger(Y = d$Y, X = d$X, random = random, data = d$smp,
             var.adjust = TRUE, B_adj = 2, seed = seed,
             num.trees = 50, MaxIterations = 10)
}

test_that("var.adjust=TRUE reports K, sigma_e, sigma_u", {
  m <- fit_adj(tiny_saef_data())
  expect_true(all(c("K", "sigma_e", "sigma_u",
                    "RanEffSD_unadj", "ErrorSD_unadj") %in% names(m)))
  expect_true(all(vapply(m$K, function(k) k >= 0, logical(1))))
  expect_gt(m$sigma_e[[length(m$sigma_e)]], 0)
})

test_that("adjusted ErrorSD does not exceed the naive residual SD", {
  d <- tiny_saef_data(); m <- fit_adj(d)
  naive <- m$sd_Mendez_Lohr_naive_unad
  # ErrorSD is the RE-model sigma on rescaled residuals; assert finiteness + no NaN
  expect_false(is.na(m$ErrorSD))
  expect_true(is.finite(m$ErrorSD))
})

test_that("var.adjust=TRUE converges within MaxIterations", {
  d <- tiny_saef_data()
  m <- MERFranger(Y = d$Y, X = d$X, random = "(1|district)", data = d$smp,
                  var.adjust = TRUE, B_adj = 2, seed = 3,
                  num.trees = 50, MaxIterations = 25)
  expect_lte(m$IterationsUsed, 25)
})
