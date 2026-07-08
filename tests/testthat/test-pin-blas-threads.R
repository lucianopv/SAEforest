test_that("pin_blas_threads single-threads a forked worker", {
  skip_on_os("windows")
  skip_on_cran()
  skip_if_not_installed("RhpcBLASctl")
  res <- parallel::mclapply(1:2, function(i) {
    pin_blas_threads()
    RhpcBLASctl::omp_get_max_threads()
  }, mc.cores = 2)
  expect_equal(unlist(res), c(1L, 1L))
})
