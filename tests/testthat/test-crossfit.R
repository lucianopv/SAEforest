# re_method = "crossfit": the random-effect step is shown predictions from a
# forest that never saw the area (simulations/PREDICTION_crossfit.md).

make_absorbed <- function(D = 12, n_i = 20, sigma_v = 0.5, n_z = 6, seed = 1) {
  set.seed(seed)
  area <- factor(rep(seq_len(D), each = n_i))
  x1 <- rnorm(D * n_i); x2 <- rnorm(D * n_i)
  v <- rnorm(D, 0, sigma_v)[as.integer(area)]
  y <- 2 * x1 + 1.5 * x2 + v + rnorm(D * n_i, 0, 1)
  # z is constant within area and carries no signal: the learner can identify the
  # area from it, which is what drives absorption.
  X <- data.frame(x1 = x1, x2 = x2)
  for (k in seq_len(n_z)) X[[paste0("z", k)]] <- runif(D)[as.integer(area)]
  list(Y = y, X = X, data = cbind(data.frame(area = area), X))
}

test_that("folds are assigned by area, never by unit", {
  grp <- factor(rep(letters[1:7], each = 3))
  f <- assign_area_folds_(grp, 3)
  expect_equal(length(f), 7L)                 # one entry per AREA
  expect_setequal(names(f), levels(grp))
  expect_true(all(f %in% 1:3))
  expect_equal(f, assign_area_folds_(grp, 3)) # deterministic
  # fewer areas than folds must not produce empty training sets
  expect_lte(max(assign_area_folds_(factor(c("a", "b")), 5)), 2L)
})

test_that("the grouping factor is read from the random formula", {
  d <- data.frame(district = factor(c("a", "b")), y = 1:2)
  expect_equal(levels(grouping_factor_(d, "(1|district)")), c("a", "b"))
  expect_error(grouping_factor_(d, "(1|missing)"), "cannot find")
  expect_error(grouping_factor_(d, "(1|a/b)"), "single grouping factor")
})

test_that("crossfit recovers a random effect that in-sample fitting absorbs", {
  d <- make_absorbed()
  common <- list(Y = d$Y, X = d$X, random = "(1|area)", data = d$data,
                 num.trees = 200, min.node.size = 5, num.threads = 1, seed = 4)

  absorbed <- do.call(MERFranger, common)
  fixed <- do.call(MERFranger, c(common, list(re_method = "crossfit", re_folds = 4)))

  # The point of the whole exercise: the default collapses, crossfit does not.
  expect_lt(absorbed$RanEffSD, 0.1)
  expect_gt(fixed$RanEffSD, 0.2)
})

test_that("crossfit leaves a well-behaved fit alone and defaults are untouched", {
  # No area-constant covariates: nothing to absorb, so both routes should agree
  # that the random effect is there.
  d <- make_absorbed(n_z = 0)
  common <- list(Y = d$Y, X = d$X, random = "(1|area)", data = d$data,
                 num.trees = 200, min.node.size = 5, num.threads = 1, seed = 4)
  a <- do.call(MERFranger, common)
  b <- do.call(MERFranger, c(common, list(re_method = "crossfit", re_folds = 4)))
  expect_gt(a$RanEffSD, 0.2)
  expect_gt(b$RanEffSD, 0.2)

  # The default argument must change nothing.
  c_explicit <- do.call(MERFranger, c(common, list(re_method = "insample")))
  expect_equal(a$RanEffSD, c_explicit$RanEffSD)
  expect_equal(a$ErrorSD, c_explicit$ErrorSD)
})

test_that("re_folds is validated", {
  d <- make_absorbed(D = 6, n_i = 10, n_z = 2)
  expect_error(
    MERFranger(Y = d$Y, X = d$X, random = "(1|area)", data = d$data,
               re_method = "crossfit", re_folds = 1, num.trees = 50, num.threads = 1),
    "re_folds")
})
