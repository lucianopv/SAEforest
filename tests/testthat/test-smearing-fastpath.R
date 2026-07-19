# Closed-form smearing indicators.
#
# The smearing set is the Minkowski sum {preds[k] + res[j]} over all k, j. All
# indicators in calc_indicat() are symmetric functions of that set, so the
# reference below (explicit materialisation) is the specification: any fast path
# must reproduce it exactly.

# Reference: materialise the full N x n smearing set and reduce it, exactly as
# point_nonLin() does today.
ref_smear <- function(preds, res, thresh, custom = NULL, select.indicator = NULL) {
  calc_indicat(as.vector(outer(preds, res, "+")),
               threshold = thresh, custom = custom,
               select.indicator = select.indicator)
}

test_that("smear_indicators matches the materialised reference for Mean", {
  preds <- c(100, 250, 400, 375.5)
  res <- c(-40, -5, 0, 12.25, 88)
  expect_equal(
    smear_indicators(preds, res, thresh = 200, select.indicator = "Mean"),
    ref_smear(preds, res, thresh = 200, select.indicator = "Mean")
  )
})

test_that("smear_indicators matches the materialised reference for Hcr", {
  preds <- c(100, 250, 400, 375.5)
  res <- c(-40, -5, 0, 12.25, 88)
  expect_equal(
    smear_indicators(preds, res, thresh = 200, select.indicator = "Hcr"),
    ref_smear(preds, res, thresh = 200, select.indicator = "Hcr")
  )
})

test_that("Hcr counts cells strictly below the threshold, excluding exact ties", {
  # every preds[k] + res[j] lands exactly on 40, 50 or 60; threshold 50 means
  # the 50-cells must NOT count (calc_indicat uses Y < threshold).
  preds <- c(10, 20, 30)
  res <- c(10, 20, 30)
  expect_equal(
    smear_indicators(preds, res, thresh = 50, select.indicator = "Hcr"),
    ref_smear(preds, res, thresh = 50, select.indicator = "Hcr")
  )
  # guard the guard: cells are {20,30,40, 30,40,50, 40,50,60}; six are strictly
  # below 50 and the two cells sitting exactly on 50 are excluded.
  expect_equal(ref_smear(preds, res, thresh = 50, select.indicator = "Hcr"),
               c(Hcr = 6 / 9))
})

test_that("smear_indicators matches the materialised reference for Pgap", {
  preds <- c(100, 250, 400, 375.5)
  res <- c(-40, -5, 0, 12.25, 88)
  expect_equal(
    smear_indicators(preds, res, thresh = 200, select.indicator = "Pgap"),
    ref_smear(preds, res, thresh = 200, select.indicator = "Pgap")
  )
})

test_that("Pgap drops negative cells from both the numerator and the denominator", {
  # calc_indicat sets y < 0 to NA and then averages with na.rm = TRUE, so
  # negative cells shrink the denominator rather than contributing zero.
  preds <- c(10, 50, 300)
  res <- c(-60, -5, 0, 40)
  expect_equal(
    smear_indicators(preds, res, thresh = 100, select.indicator = "Pgap"),
    ref_smear(preds, res, thresh = 100, select.indicator = "Pgap")
  )
  # a zero-contribution treatment would divide by 12 instead of 11
  expect_false(isTRUE(all.equal(
    unname(ref_smear(preds, res, thresh = 100, select.indicator = "Pgap")),
    sum(pmax(0, 100 - pmax(0, as.vector(outer(preds, res, "+")))) / 100) / 12
  )))
})

test_that("Pgap is NaN when every smeared cell is negative", {
  preds <- c(-100, -200)
  res <- c(-10, -20, -30)
  expect_equal(
    smear_indicators(preds, res, thresh = 50, select.indicator = "Pgap"),
    ref_smear(preds, res, thresh = 50, select.indicator = "Pgap")
  )
  expect_true(is.nan(unname(ref_smear(preds, res, thresh = 50, select.indicator = "Pgap"))))
})

# --- fallback routing --------------------------------------------------------
# Only Mean/Hcr/Pgap have closed forms. Anything requiring the sorted smearing
# set (Gini, quantiles, Qsr) or an arbitrary user function must still be served
# by materialising, bit-identically to the current implementation.

test_that("smear_indicators falls back to materialisation when no selection is given", {
  preds <- c(100, 250, 400, 375.5)
  res <- c(-40, -5, 0, 12.25, 88)
  expect_equal(
    smear_indicators(preds, res, thresh = 200),
    ref_smear(preds, res, thresh = 200)
  )
})

test_that("smear_indicators falls back when a sort-requiring indicator is selected", {
  preds <- c(100, 250, 400, 375.5)
  res <- c(-40, -5, 0, 12.25, 88)
  for (sel in list("Gini", "Median", "Quant10", "Qsr", c("Mean", "Gini"))) {
    expect_equal(
      smear_indicators(preds, res, thresh = 200, select.indicator = sel),
      ref_smear(preds, res, thresh = 200, select.indicator = sel),
      info = paste(sel, collapse = "+")
    )
  }
})

test_that("smear_indicators falls back when custom indicators are supplied", {
  preds <- c(100, 250, 400, 375.5)
  res <- c(-40, -5, 0, 12.25, 88)
  custom <- list(my_p10 = function(Y, threshold) mean(Y[Y <= quantile(Y, 0.1)]))
  # select the custom column itself, otherwise it is computed and discarded and
  # the assertion would pass without the fallback ever being needed
  expect_equal(
    smear_indicators(preds, res, thresh = 200, custom = custom,
                     select.indicator = c("Mean", "my_p10")),
    ref_smear(preds, res, thresh = 200, custom = custom,
              select.indicator = c("Mean", "my_p10"))
  )
  expect_equal(
    smear_indicators(preds, res, thresh = 200, custom = custom),
    ref_smear(preds, res, thresh = 200, custom = custom)
  )
})

# --- end-to-end --------------------------------------------------------------

test_that("SAEforest_model indicators do not depend on which ones are selected", {
  skip_on_cran()
  d <- tiny_saef_data()
  run <- function(sel) {
    set.seed(2026)
    SAEforest_model(Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp,
                    pop_data = d$pop, meanOnly = FALSE, MSE = "none",
                    num.trees = 50, mtry = 3, select.indicator = sel)$Indicators
  }
  full <- run(NULL)
  sub <- run(c("Mean", "Hcr", "Pgap"))
  expect_equal(sub, full[, c("district", "Mean", "Hcr", "Pgap")])
})
