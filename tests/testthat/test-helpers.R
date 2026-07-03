test_that("adjust_ErrorSD_ returns a finite non-negative K", {
  d <- tiny_saef_data()
  set.seed(5)
  rf <- ranger::ranger(x = d$X, y = d$Y, num.trees = 50)
  K <- adjust_ErrorSD_(Y = d$Y, X = d$X, smp_data = d$smp, rf = rf,
                       B = 2, num.trees = 50)
  expect_length(K, 1)
  expect_true(is.finite(K))
  expect_gte(K, 0)
})

test_that("adjust_ErrorSD clamps to 0 (not NaN) and warns when K exceeds ErrorSD^2", {
  d <- tiny_saef_data()
  fake <- list(Forest = ranger::ranger(x = d$X, y = d$Y, num.trees = 50, seed = 5),
               ErrorSD = 1e-6)  # tiny SD forces the clamp
  expect_warning(
    out <- adjust_ErrorSD(Y = d$Y, X = d$X, smp_data = d$smp, mod = fake,
                          B = 2, num.trees = 50),
    "exceeds"
  )
  expect_false(is.nan(out)); expect_gte(out, 0)
})

test_that("calc_indicat computes a correct head count ratio", {
  y <- c(1, 2, 3, 4, 100)
  ind <- calc_indicat(y, threshold = 3, custom = NULL, select.indicator = "Hcr")
  expect_equal(unname(ind[["Hcr"]]), mean(y < 3))
})

test_that("ran_comp does not return NaN for a constant-residual domain", {
  d <- tiny_saef_data()
  mod <- MERFranger(Y = d$Y, X = d$X, random = "(1|district)", data = d$smp,
                    seed = 5, num.trees = 50)
  rc <- ran_comp(mod = mod, smp_data = d$smp, Y = d$Y,
                 dName = "district", ADJsd = mod$ErrorSD)
  expect_false(any(is.nan(rc$forest_res)))
  expect_false(any(is.nan(rc$ran_effs)))
})
