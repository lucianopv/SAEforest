# force_serial_threads() gained `worker.threads`, letting each forked MSE worker
# run multi-threaded ranger instead of the pinned single thread. Memory scales
# with the number of WORKERS (private allocations) while threads share one
# address space, so trading processes for threads at a fixed CPU count cuts
# memory without losing parallelism. Default stays 1: opting in is deliberate,
# because cores x threads over the allocation oversubscribes the node, which is
# why this function pinned to 1 in the first place.

test_that("default is unchanged: workers still pinned to one thread", {
  d <- force_serial_threads(4L, num.trees = 50)
  expect_identical(d$num.threads, 1L)
  expect_identical(d$num.trees, 50)
})

test_that("cores == 1 is untouched (no fork, nothing to pin)", {
  d <- force_serial_threads(1L, num.trees = 50)
  expect_null(d$num.threads)
})

test_that("worker.threads sets ranger's thread count inside workers", {
  d <- force_serial_threads(2L, num.trees = 50, worker.threads = 4L)
  expect_identical(d$num.threads, 4L)
  expect_false("worker.threads" %in% names(d))   # must not leak through to ranger
})

test_that("cores x worker.threads beyond the allocation errors loudly", {
  withr::with_envvar(c(SLURM_CPUS_PER_TASK = "4"), {
    expect_error(force_serial_threads(2L, worker.threads = 4L), "oversubscribe")
    expect_silent(force_serial_threads(2L, worker.threads = 2L))
  })
})

test_that("no allocation set means no guard, but still applies the count", {
  withr::with_envvar(c(SLURM_CPUS_PER_TASK = NA), {
    d <- force_serial_threads(2L, worker.threads = 8L)
    expect_identical(d$num.threads, 8L)
  })
})

test_that("worker.threads must be a positive integer", {
  expect_error(force_serial_threads(2L, worker.threads = 0L), "positive integer")
})
