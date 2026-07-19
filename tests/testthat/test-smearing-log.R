# Log-path (multiplicative) smearing. Cells are z[k,j] = exp(preds[k]) * exp(res[j]).
# point_nonLin() maps non-finite cells to NA before reducing, so the na.rm
# denominators are the count of FINITE cells -- the reference below reproduces
# that exactly and is the specification.

ref_log <- function(preds, res, thresh, custom = NULL, select.indicator = NULL) {
  v <- as.vector(outer(exp(preds), exp(res), "*"))
  v[!is.finite(v)] <- NA
  calc_indicat(v, threshold = thresh, custom = custom,
               select.indicator = select.indicator)
}

lp <- c(9.1, 10.4, 8.25, 11.0, 9.75)     # log-scale predictions
lr <- c(-0.8, -0.1, 0, 0.35, 1.2, 0.05)  # log-scale residuals
TH <- 20000

test_that("log-path Mean matches the materialised reference", {
  expect_equal(
    smear_indicators(lp, lr, thresh = TH, select.indicator = "Mean",
                     transformation = "log"),
    ref_log(lp, lr, thresh = TH, select.indicator = "Mean")
  )
})

test_that("log-path Hcr matches the materialised reference", {
  expect_equal(
    smear_indicators(lp, lr, thresh = TH, select.indicator = "Hcr",
                     transformation = "log"),
    ref_log(lp, lr, thresh = TH, select.indicator = "Hcr")
  )
})

test_that("log-path Pgap matches the materialised reference", {
  expect_equal(
    smear_indicators(lp, lr, thresh = TH, select.indicator = "Pgap",
                     transformation = "log"),
    ref_log(lp, lr, thresh = TH, select.indicator = "Pgap")
  )
})

test_that("log-path handles all three together and across thresholds", {
  sel <- c("Mean", "Hcr", "Pgap")
  for (th in c(1, 5000, 20000, 1e6)) {
    expect_equal(
      smear_indicators(lp, lr, thresh = th, select.indicator = sel,
                       transformation = "log"),
      ref_log(lp, lr, thresh = th, select.indicator = sel),
      info = paste("threshold", th)
    )
  }
})

test_that("log-path falls back to materialisation when cells overflow to non-finite", {
  # exp(750) is Inf; the reference turns those cells into NA, which changes the
  # na.rm denominators. The closed forms assume every cell is finite, so this
  # must route to the materialised path rather than compute a wrong denominator.
  big <- c(9.1, 750, 10.0)
  expect_true(any(!is.finite(as.vector(outer(exp(big), exp(lr), "*")))))
  sel <- c("Mean", "Hcr", "Pgap")
  expect_equal(
    smear_indicators(big, lr, thresh = TH, select.indicator = sel,
                     transformation = "log"),
    ref_log(big, lr, thresh = TH, select.indicator = sel)
  )
})

test_that("log-path falls back for Gini, Qsr and custom indicators", {
  for (sel in list("Gini", "Qsr", c("Mean", "Gini"))) {
    expect_equal(
      smear_indicators(lp, lr, thresh = TH, select.indicator = sel,
                       transformation = "log"),
      ref_log(lp, lr, thresh = TH, select.indicator = sel),
      info = paste(sel, collapse = "+")
    )
  }
  custom <- list(cv = function(Y, threshold) sd(Y, na.rm = TRUE) / mean(Y, na.rm = TRUE))
  expect_equal(
    smear_indicators(lp, lr, thresh = TH, custom = custom, transformation = "log"),
    ref_log(lp, lr, thresh = TH, custom = custom)
  )
})
