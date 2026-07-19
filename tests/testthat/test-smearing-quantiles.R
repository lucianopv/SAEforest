# r-th order statistic of the Minkowski sum {A[k] + B[j]}, without materialising
# it. With A and B sorted the implied matrix is sorted along both axes, so the
# order statistic can be found by bisection on value plus a counting pass.

test_that("minkowski_order_stat returns the r-th smallest of the Minkowski sum", {
  A <- c(100, 250, 400, 375.5)
  B <- c(-40, -5, 0, 12.25, 88)
  full <- sort(as.vector(outer(A, B, "+")))
  As <- sort(A); Bs <- sort(B)
  for (r in seq_along(full)) {
    expect_equal(minkowski_order_stat(As, Bs, r), full[r], info = paste("r =", r))
  }
})

test_that("minkowski_order_stat handles heavy ties", {
  # many cells collapse onto identical values
  A <- c(0, 10, 20, 30)
  B <- c(0, 10, 20, 30)
  full <- sort(as.vector(outer(A, B, "+")))
  As <- sort(A); Bs <- sort(B)
  for (r in seq_along(full)) {
    expect_equal(minkowski_order_stat(As, Bs, r), full[r], info = paste("r =", r))
  }
})

test_that("minkowski_order_stat handles negative and duplicated inputs", {
  set.seed(11)
  A <- c(rnorm(6, -5, 3), 0, 0)
  B <- c(rnorm(5, 2, 4), 0)
  full <- sort(as.vector(outer(A, B, "+")))
  As <- sort(A); Bs <- sort(B)
  for (r in c(1L, 5L, 13L, 24L, length(full))) {
    expect_equal(minkowski_order_stat(As, Bs, r), full[r], info = paste("r =", r))
  }
})

# --- integration -------------------------------------------------------------

ref_smear2 <- function(preds, res, thresh, custom = NULL, select.indicator = NULL) {
  calc_indicat(as.vector(outer(preds, res, "+")), threshold = thresh,
               custom = custom, select.indicator = select.indicator)
}

test_that("smear_indicators matches the reference for every quantile", {
  preds <- c(100, 250, 400, 375.5, 12, 908.25)
  res <- c(-40, -5, 0, 12.25, 88, 3.5, 77)
  for (q in c("Quant10", "Quant25", "Median", "Quant75", "Quant90")) {
    expect_equal(
      smear_indicators(preds, res, thresh = 200, select.indicator = q),
      ref_smear2(preds, res, thresh = 200, select.indicator = q), info = q
    )
  }
  sel <- c("Mean", "Median", "Hcr", "Quant90")
  expect_equal(
    smear_indicators(preds, res, thresh = 200, select.indicator = sel),
    ref_smear2(preds, res, thresh = 200, select.indicator = sel)
  )
})

test_that("Gini and Qsr still fall back to the materialised path", {
  preds <- c(100, 250, 400, 375.5, 12, 908.25)
  res <- c(-40, -5, 0, 12.25, 88, 3.5, 77)
  for (sel in list("Gini", "Qsr", c("Median", "Gini"), c("Quant10", "Qsr"))) {
    expect_equal(
      smear_indicators(preds, res, thresh = 200, select.indicator = sel),
      ref_smear2(preds, res, thresh = 200, select.indicator = sel),
      info = paste(sel, collapse = "+")
    )
  }
})

test_that("quantiles agree with the reference on the selection path (above the size guard)", {
  # smaller sets route to materialisation; this one must exceed
  # smear_select_min_cells so the order-statistic path is the one under test
  set.seed(4)
  preds <- rnorm(500, 20000, 6000)
  res <- rnorm(200, 0, 4000)
  expect_gte(length(preds) * length(res), smear_select_min_cells)
  sel <- c("Mean", "Hcr", "Pgap", "Quant10", "Quant25", "Median", "Quant75", "Quant90")
  expect_equal(
    smear_indicators(preds, res, thresh = 10000, select.indicator = sel),
    ref_smear2(preds, res, thresh = 10000, select.indicator = sel)
  )
})
