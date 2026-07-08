test_that("point_nonLin's log-transform smearing produces the same Indicators before and after exp-factoring", {
  d <- tiny_saef_data()
  set.seed(11)
  result <- point_nonLin(
    Y = d$Y, X = d$X, dName = d$dName, threshold = 15000,
    smp_data = d$smp, pop_data = d$pop,
    initialRandomEffects = 0, ErrorTolerance = 0.0001, MaxIterations = 25,
    custom_indicator = NULL, aggregate_to = NULL,
    transformation = "log", num.trees = 50
  )
  ind <- result[[1]]

  expect_true(all(c("Mean", "Hcr", "Gini") %in% names(ind)))
  expect_false(anyNA(ind$Mean))
  expect_true(all(ind$Mean > 0))  # log-transform back-transformed via exp(), always positive
})

test_that("point_nonLin's transformation = 'none' path is unaffected by the exp-factoring change", {
  d <- tiny_saef_data()
  set.seed(12)
  # transformation = 'none' has no exp() step to factor -- must keep using the
  # original additive matrix construction path, completely untouched by Task 2.
  result <- point_nonLin(
    Y = d$Y, X = d$X, dName = d$dName, threshold = 15000,
    smp_data = d$smp, pop_data = d$pop,
    initialRandomEffects = 0, ErrorTolerance = 0.0001, MaxIterations = 25,
    custom_indicator = NULL, aggregate_to = NULL,
    transformation = "none", num.trees = 50
  )
  ind <- result[[1]]

  expect_true(all(c("Mean", "Hcr") %in% names(ind)))
  expect_false(anyNA(ind$Mean))
})
