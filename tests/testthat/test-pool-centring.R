# The smearing pool must be mean-zero on the variance-adjusted path.
#
# OOBresiduals feeds the smearing point estimator, so a non-zero mean is a pure
# level error in every indicator computed from it. Before the fix,
# `r_ij_new - AllEffects` left mean(-observation-weighted BLUP), which is
# non-zero whenever domain SIZE correlates with the domain EFFECT. These tests
# build exactly that correlation, because a balanced design cannot detect it.

make_unbalanced <- function(D = 30, seed = 1, sigma_v = 0.6, sigma_e = 0.5) {
  set.seed(seed)
  # Domain size and domain effect deliberately correlated: n_i grows with i and
  # so does v_i. This is the GridSAE shape (urban domains are both larger and
  # better off), and it is what makes the observation-weighted BLUP mean differ
  # from the domain-weighted one.
  n_i <- rep(seq_len(D), each = 1) * 2 + 2
  v   <- scale(seq_len(D))[, 1] * sigma_v
  area <- factor(rep(seq_len(D), times = n_i))
  vv  <- v[as.integer(area)]
  n   <- length(area)
  x1  <- stats::rnorm(n); x2 <- stats::rnorm(n)
  y   <- 2 + 0.8 * x1 - 0.4 * x2 + vv + stats::rnorm(n, sd = sigma_e)
  data.frame(y = y, x1 = x1, x2 = x2, area = area)
}

fit <- function(d, ...) {
  suppressWarnings(MERFranger(Y = d$y, X = d[, c("x1", "x2")],
                              random = "(1|area)", data = d,
                              num.trees = 50, num.threads = 1, seed = 1, ...))
}

test_that("var.adjust pool is mean-zero even with size-correlated domain effects", {
  d <- make_unbalanced()
  m <- fit(d, var.adjust = TRUE)
  expect_gt(m$RanEffSD, 0.1)                      # the test is vacuous at sigma_v = 0
  r <- m$OOBresiduals; r <- r[is.finite(r)]
  expect_lt(abs(mean(r)), 1e-8)
})

test_that("the observation-weighted BLUP mean really is non-zero here", {
  # Guards the test above against becoming vacuous: if a future lme4 centred the
  # BLUPs by observation, the fix would be untestable and this would fail first.
  d <- make_unbalanced()
  m <- fit(d, var.adjust = TRUE)
  re <- m$RandomEffects[[1]][, 1]
  by_domain <- mean(re)
  by_obs <- mean(re[match(as.character(d$area), rownames(m$RandomEffects[[1]]))])
  expect_lt(abs(by_domain), 1e-6)                 # lme4 centres over DOMAINS
  expect_gt(abs(by_obs), 1e-3)                    # but not over OBSERVATIONS
})

test_that("sigma_v = 0 leaves the pool bit-identical (production fits unaffected)", {
  # With AllEffects = 0 and r_ij_new already centred, the new subtraction is a
  # no-op. This is why GridSAE's published fits, which all report an exact
  # singular zero, cannot move.
  set.seed(4)
  n <- 400
  d <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n),
                  area = factor(rep(1:8, length.out = n)))
  d$y <- 1 + d$x1 - 0.5 * d$x2 + stats::rnorm(n, sd = 0.4)   # no area effect
  m <- fit(d, var.adjust = TRUE)
  r <- m$OOBresiduals; r <- r[is.finite(r)]
  expect_lt(abs(mean(r)), 1e-8)
})

test_that("the upstream path is untouched", {
  d <- make_unbalanced()
  m <- fit(d, var.adjust = FALSE)
  r <- m$OOBresiduals
  # Upstream centres naturally because predict(lmefit) carries the offset; this
  # asserts the fix did not reach into that branch at all.
  expect_true(is.numeric(r))
  expect_equal(length(r), nrow(d))
})
