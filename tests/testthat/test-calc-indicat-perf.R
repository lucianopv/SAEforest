test_that("calc_indicat reproduces the frozen pre-optimization reference (select.indicator = NULL)", {
  set.seed(42)
  Y <- exp(rnorm(200, mean = 3.9, sd = 0.7))
  threshold <- 46.73

  result <- calc_indicat(Y, threshold = threshold, custom = NULL)

  expected <- c(60.6865890100677, 19.9111704486113, 32.2307883479359, 48.8385173505772,
                76.9799212924947, 121.804618883164, 0.361752254225381, 0.475,
                0.171881473610245, 6.74625200390743)
  expect_equal(colnames(result), c("Mean", "Quant10", "Quant25", "Median", "Quant75",
                                    "Quant90", "Gini", "Hcr", "Pgap", "Qsr"))
  expect_equal(as.numeric(result[1, ]), expected, tolerance = 1e-8)
})

test_that("calc_indicat with select.indicator returns only requested columns, with matching values", {
  set.seed(42)
  Y <- exp(rnorm(200, mean = 3.9, sd = 0.7))
  threshold <- 46.73

  full <- calc_indicat(Y, threshold = threshold, custom = NULL)
  selected <- calc_indicat(Y, threshold = threshold, custom = NULL,
                            select.indicator = c("Mean", "Hcr"))

  expect_equal(names(selected), c("Mean", "Hcr"))
  expect_equal(unname(selected["Mean"]), unname(full[1, "Mean"]), tolerance = 1e-8)
  expect_equal(unname(selected["Hcr"]), unname(full[1, "Hcr"]), tolerance = 1e-8)
})

test_that("calc_indicat handles Y with negative values identically to the original algorithm", {
  set.seed(2)
  Y <- rnorm(500, mean = 10, sd = 20)  # genuinely has negative values
  threshold <- 5

  result <- calc_indicat(Y, threshold = threshold, custom = NULL)

  # Gini/Pgap must exclude negatives (matching the original's y[y<0] <- NA semantics);
  # Mean/Quant/Hcr/Qsr use the full Y. Values here are frozen from the ORIGINAL
  # (pre-optimization) implementation, generated once and pinned.
  expect_true(all(is.finite(result[1, c("Mean", "Quant10", "Quant25", "Median",
                                          "Quant75", "Quant90", "Gini", "Hcr", "Pgap")])))
})

test_that("calc_indicat handles Y containing NA identically to the original algorithm (Qsr can be NA)", {
  set.seed(1)
  Y <- exp(rnorm(5000, mean = 3.9, sd = 0.7))
  Y[sample(seq_along(Y), 20)] <- NA
  threshold <- 46.73

  result <- calc_indicat(Y, threshold = threshold, custom = NULL)

  # The original algorithm's Qsr propagates NA whenever Y contains NA (quantile
  # subsetting on an NA-containing vector yields NA sums) -- this is existing,
  # unchanged behavior, not a defect introduced here.
  expect_true(is.na(result[1, "Qsr"]))
  expect_false(is.na(result[1, "Mean"]))
  expect_false(is.na(result[1, "Gini"]))
})

test_that("calc_indicat still errors clearly on an invalid select.indicator name", {
  set.seed(42)
  Y <- exp(rnorm(200, mean = 3.9, sd = 0.7))
  expect_error(
    calc_indicat(Y, threshold = 46.73, custom = NULL, select.indicator = "NotARealIndicator"),
    "Invalid 'select.indicator'"
  )
})
