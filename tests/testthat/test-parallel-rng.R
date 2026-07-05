test_that("cores = 1 is serial and does not touch the RNG kind", {
  before <- RNGkind()
  g <- function(cl) pbapply::pbsapply(1:4, function(i) i^2, cl = cl)
  a <- with_parallel_rng(1, g)
  b <- pbapply::pbsapply(1:4, function(i) i^2, cl = NULL)
  expect_equal(a, b)
  expect_identical(RNGkind(), before)
})

test_that("cores > 1 is reproducible under a fixed outer seed and restores the RNG", {
  skip_on_cran(); skip_on_os("windows")
  before <- RNGkind()
  g <- function(cl) pbapply::pbsapply(1:8, function(i) stats::rnorm(1), cl = cl)
  set.seed(10); a <- with_parallel_rng(2, g)
  set.seed(10); b <- with_parallel_rng(2, g)
  expect_equal(a, b)                    # same (seed, cores) -> same result
  expect_identical(RNGkind(), before)   # RNG kind restored after
})

test_that("force_serial_threads injects num.threads=1 only when parallel", {
  expect_equal(force_serial_threads(1, num.trees = 50), list(num.trees = 50))
  expect_equal(force_serial_threads(2, num.trees = 50)$num.threads, 1L)
  expect_warning(force_serial_threads(2, num.threads = 4), "num.threads")
  expect_silent(force_serial_threads(2, num.trees = 50))            # no conflict -> no warning
})
