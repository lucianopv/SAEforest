# Integer-overflow guard: length(preds) * length(res) exceeds .Machine$integer.max
# for a large domain (GridSAE's biggest district is 280,632 units x 15,270
# residuals = 4.3e9 cells). Integer arithmetic there yields NA and silently
# poisons Hcr/Pgap.
#
# Replicating every prediction k times leaves all indicators unchanged -- the
# smearing set's value multiset is scaled uniformly -- so a small case with a
# known answer can be blown up past the integer limit and must give the SAME
# numbers. Cheap to run precisely because the fast path never materialises.

test_that("Hcr/Pgap survive a domain whose cell count exceeds integer.max (log path)", {
  set.seed(3)
  p0 <- rnorm(12, 9.8, 0.6)
  res <- rnorm(15000, 0, 0.45)
  small <- smear_indicators(p0, res, thresh = 18000,
                            select.indicator = c("Mean", "Hcr", "Pgap"),
                            transformation = "log")

  p_big <- rep(p0, 20000)                       # 240,000 x 15,000 = 3.6e9 cells
  expect_gt(length(p_big) * as.numeric(length(res)), .Machine$integer.max)
  big <- smear_indicators(p_big, res, thresh = 18000,
                          select.indicator = c("Mean", "Hcr", "Pgap"),
                          transformation = "log")

  expect_true(all(is.finite(big)))
  expect_equal(big, small)
})

test_that("Hcr/Pgap survive a domain whose cell count exceeds integer.max (additive path)", {
  set.seed(3)
  p0 <- rnorm(12, 20000, 6000)
  res <- rnorm(15000, 0, 4000)
  small <- smear_indicators(p0, res, thresh = 10000,
                            select.indicator = c("Mean", "Hcr", "Pgap"))

  p_big <- rep(p0, 20000)
  expect_gt(length(p_big) * as.numeric(length(res)), .Machine$integer.max)
  big <- smear_indicators(p_big, res, thresh = 10000,
                          select.indicator = c("Mean", "Hcr", "Pgap"))

  expect_true(all(is.finite(big)))
  expect_equal(big, small)
})

test_that("quantile selection survives a cell count exceeding integer.max", {
  set.seed(3)
  p0 <- rnorm(12, 20000, 6000)
  res <- rnorm(15000, 0, 4000)
  p_big <- rep(p0, 20000)
  expect_gt(length(p_big) * as.numeric(length(res)), .Machine$integer.max)
  big <- smear_indicators(p_big, res, thresh = 10000,
                          select.indicator = c("Median", "Quant90"))
  expect_true(all(is.finite(big)))
})
