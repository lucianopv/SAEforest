# Fixed-weight variance bias correction (Mendez & Lohr 2011, Procedure 2),
# opt-in via MERFranger(var.adjust = TRUE, adjust_method = "weights").

# make_toy() lives in helper-fixtures.R (shared with test-MERFranger.R).

test_that("the weight matrix reproduces ranger's own OOB predictions", {
  d <- make_toy()
  rf <- ranger::ranger(x = d$X, y = d$Y, num.trees = 200, min.node.size = 5,
                       keep.inbag = TRUE, seed = 7, num.threads = 1)
  W <- oob_weight_matrix_(rf, d$X)
  ok <- attr(W, "n_oob") > 0 & is.finite(rf$predictions)

  # This is the property the whole procedure rests on: the forest IS this matrix.
  expect_equal(as.vector(W %*% d$Y)[ok], rf$predictions[ok], tolerance = 1e-10)
  expect_equal(rowSums(W)[ok], rep(1, sum(ok)), tolerance = 1e-10)
  expect_true(all(diag(W) == 0))          # never predicts itself
  expect_true(all(W >= 0))
})

test_that("a forest without keep.inbag is refused rather than guessed at", {
  d <- make_toy()
  rf <- ranger::ranger(x = d$X, y = d$Y, num.trees = 50, seed = 7, num.threads = 1)
  expect_error(oob_weight_matrix_(rf, d$X), "keep.inbag")
})

test_that("K from fixed weights is finite, positive and below the refit K", {
  d <- make_toy()
  rf <- ranger::ranger(x = d$X, y = d$Y, num.trees = 200, min.node.size = 5,
                       keep.inbag = TRUE, seed = 7, num.threads = 1)
  set.seed(11)
  kw <- adjust_ErrorSD_weights_(Y = d$Y, X = d$X, rf = rf, B = 30)
  set.seed(11)
  kr <- adjust_ErrorSD_(Y = d$Y, X = d$X, smp_data = d$data, rf = rf, B = 30,
                        adj_tol = 0, num.trees = 200, min.node.size = 5,
                        num.threads = 1)

  expect_true(is.finite(kw$K) && kw$K > 0)
  expect_identical(kw$m, 30L)
  # The paper's point: refitting a forest per draw overcorrects relative to
  # holding the weights fixed.
  expect_lt(kw$K, kr$K)
})

test_that("MERFranger runs the weights route and leaves the default alone", {
  d <- make_toy()
  common <- list(Y = d$Y, X = d$X, random = "(1|area)", data = d$data,
                 num.trees = 100, min.node.size = 5, num.threads = 1, seed = 3)

  m_ref <- do.call(MERFranger, c(common, list(var.adjust = TRUE, B_adj = 10, adj_tol = 0)))
  m_wts <- do.call(MERFranger, c(common, list(var.adjust = TRUE, B_adj = 10,
                                              adjust_method = "weights")))

  expect_s3_class(m_wts, "MERFranger")
  expect_true(all(is.finite(unlist(m_wts$K))))
  expect_true(is.finite(m_wts$ErrorSD) && m_wts$ErrorSD > 0)
  # Same estimand, different route: the answers should differ but stay comparable.
  expect_false(isTRUE(all.equal(unlist(m_ref$K), unlist(m_wts$K))))
  expect_equal(m_wts$ErrorSD, m_ref$ErrorSD, tolerance = 0.5)

  # keep.inbag is switched on for this route, and only this route.
  expect_false(is.null(m_wts$Forest$inbag.counts))
  expect_true(is.null(m_ref$Forest$inbag.counts))
})

test_that("the weights route is refused where it would be meaningless", {
  d <- make_toy()
  common <- list(Y = d$Y, X = d$X, random = "(1|area)", data = d$data,
                 num.trees = 50, num.threads = 1)

  # Ignored, not rejected, when there is no correction to compute: wrappers such
  # as point_mean() forward their dots to MERFranger without var.adjust.
  m <- do.call(MERFranger, c(common, list(var.adjust = FALSE, adjust_method = "weights")))
  expect_s3_class(m, "MERFranger")
  expect_true(is.null(m$Forest$inbag.counts))   # the refit path, so no inbag kept
  # K_fixed means K is not estimated here at all, so adjust_method is moot and
  # ignored -- adjust_mse = "plugin" forwards both.
  m_plug <- do.call(MERFranger, c(common, list(var.adjust = TRUE, B_adj = 5,
                                               adjust_method = "weights", K_fixed = 0.1)))
  expect_s3_class(m_plug, "MERFranger")
  expect_equal(unlist(m_plug$K)[[1]], 0.1)
  expect_error(
    do.call(MERFranger, c(common, list(var.adjust = TRUE, adjust_method = "weights",
                                       keep.inbag = FALSE))),
    "keep.inbag = TRUE")
})

test_that("var.adjust = FALSE is untouched by the new argument", {
  d <- make_toy()
  common <- list(Y = d$Y, X = d$X, random = "(1|area)", data = d$data,
                 num.trees = 100, min.node.size = 5, num.threads = 1, seed = 5)
  a <- do.call(MERFranger, c(common, list(var.adjust = FALSE)))
  b <- do.call(MERFranger, c(common, list(var.adjust = FALSE, adjust_method = "refit")))
  expect_equal(a$ErrorSD, b$ErrorSD)
  expect_equal(a$RanEffSD, b$RanEffSD)
  expect_equal(unname(predict(a$Forest, d$X)$predictions),
               unname(predict(b$Forest, d$X)$predictions))
})

test_that("matrix-free OOB predictions equal the weight matrix times the outcomes", {
  d <- make_toy(); set.seed(1)
  rf <- ranger::ranger(x = d$X, y = d$Y, num.trees = 30, keep.inbag = TRUE, num.threads = 1)
  Ys <- matrix(stats::rnorm(nrow(d$X) * 4), ncol = 4)
  W <- oob_weight_matrix_(rf, d$X); ok <- attr(W, "n_oob") > 0
  expect_equal(oob_predict_boot_(rf, d$X, Ys)[ok, ], (W %*% Ys)[ok, ], tolerance = 1e-12)
})
test_that("K is the dense route's K, draw for draw", {
  d <- make_toy(); set.seed(1)
  rf <- ranger::ranger(x = d$X, y = d$Y, num.trees = 30, keep.inbag = TRUE, num.threads = 1)
  set.seed(9); k_new <- adjust_ErrorSD_weights_(d$Y, d$X, rf, B = 7)$K
  set.seed(9); k_ref <- dense_K_reference(d$Y, d$X, rf, B = 7)   # today's body, kept in helper-fixtures.R
  expect_equal(k_new, k_ref$K, tolerance = 1e-12)
})
test_that("memory stays bounded at n = 12,000", {
  set.seed(2); n <- 12000; X <- data.frame(a = stats::rnorm(n), b = stats::rnorm(n)); y <- X$a + stats::rnorm(n)
  rf <- ranger::ranger(x = X, y = y, num.trees = 10, keep.inbag = TRUE)
  invisible(gc(reset = TRUE)); adjust_ErrorSD_weights_(y, X, rf, B = 5); g <- gc()
  expect_lt(g[2, 6], 500)   # Vcells "max used" (Mb); the dense W alone is 1,099 Mb
})
