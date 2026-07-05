test_that("adjust_ErrorSD_ returns a finite non-negative K", {
  d <- tiny_saef_data()
  set.seed(5)
  rf <- ranger::ranger(x = d$X, y = d$Y, num.trees = 50)
  res <- adjust_ErrorSD_(Y = d$Y, X = d$X, smp_data = d$smp, rf = rf,
                       B = 2, num.trees = 50)
  expect_named(res, c("K", "m"))
  expect_length(res$K, 1)
  expect_true(is.finite(res$K))
  expect_gte(res$K, 0)
  expect_equal(res$m, 2)
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

test_that("calc_indicat rejects an unknown select.indicator with a clear error", {
  y <- c(1, 2, 3, 4, 100)
  err <- tryCatch(
    calc_indicat(y, threshold = 3, custom = NULL, select.indicator = "Nope"),
    error = function(e) conditionMessage(e)
  )
  # A helpful message that names the offender and lists valid indicators,
  # NOT the cryptic base-R "subscript out of bounds".
  expect_match(err, "select.indicator")
  expect_match(err, "Hcr")
  expect_false(grepl("subscript out of bounds", err))
})

test_that("ran_comp stays finite when within-domain residual SD is zero (single-obs domains)", {
  # Single-observation domains force forest_eij (the within-domain residual
  # deviation) to be exactly 0 for every row, so sd(forest_eij) == 0. The
  # old, unguarded ran_comp computed (0 / 0) * ADJsd == NaN in this case.
  set.seed(1)
  smp <- data.frame(district = factor(c("A", "B", "C")),
                    x1 = rnorm(3), y = rnorm(3))
  rf <- ranger::ranger(y = smp$y, x = smp[, "x1", drop = FALSE], num.trees = 10)
  mod <- list(Forest = rf, RanEffSD = 1)
  rc <- ran_comp(mod = mod, smp_data = smp, Y = smp$y,
                 dName = "district", ADJsd = 1)
  expect_false(any(is.nan(rc$forest_res)))
  expect_false(any(is.nan(rc$ran_effs)))
  expect_true(all(is.finite(rc$forest_res)))
})
