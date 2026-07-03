test_that("same seed gives identical point estimates", {
  d <- tiny_saef_data()
  run <- function() SAEforest_model(
    Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp, pop_data = d$pop,
    meanOnly = FALSE, MSE = "none", num.trees = 50, mtry = 3, seed = 99
  )$Indicators
  expect_equal(run(), run())
})

test_that("different seeds change the estimates", {
  d <- tiny_saef_data()
  one <- SAEforest_model(Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp,
    pop_data = d$pop, meanOnly = FALSE, MSE = "none",
    num.trees = 50, mtry = 3, seed = 1)$Indicators
  two <- SAEforest_model(Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp,
    pop_data = d$pop, meanOnly = FALSE, MSE = "none",
    num.trees = 50, mtry = 3, seed = 2)$Indicators
  expect_false(isTRUE(all.equal(one, two)))
})
