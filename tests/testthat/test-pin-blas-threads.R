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

test_that("cores > 1 nonLin MSE bootstrap completes and returns MSE estimates", {
  skip_on_os("windows")
  skip_on_cran()
  skip_if_not_installed("emdi")
  skip_if_not_installed("RhpcBLASctl")
  d <- tiny_saef_data()
  out <- SAEforest_model(
    Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp, pop_data = d$pop,
    MSE = "nonparametric", B = 2, num.trees = 50, threshold = median(d$Y),
    seed = 1, cores = 2
  )
  expect_s3_class(out, "SAEforest")
  expect_false(is.null(out$MSE_Estimates))
})
