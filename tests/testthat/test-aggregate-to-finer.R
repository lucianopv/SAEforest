test_that("aggregate_to rejects a column that does not nest within dName", {
  d <- tiny_saef_data_with_subunit()
  d$pop$bad_group <- sample(1:3, nrow(d$pop), replace = TRUE)  # independent of district

  expect_error(
    SAEforest_model(
      Y = d$Y, X = d$X, dName = d$dName,
      smp_data = d$smp, pop_data = d$pop,
      meanOnly = FALSE, MSE = "none", num.trees = 50, mtry = 2,
      aggregate_to = "bad_group"
    ),
    "does not nest within"
  )
})

test_that("aggregate_to does not change fitted model parameters (fixef, VarCorr, RanEffSD)", {
  d <- tiny_saef_data_with_subunit()

  set.seed(1)
  plain <- SAEforest_model(
    Y = d$Y, X = d$X, dName = d$dName,
    smp_data = d$smp, pop_data = d$pop,
    meanOnly = FALSE, MSE = "none", num.trees = 50, mtry = 2
  )

  set.seed(1)
  finer <- SAEforest_model(
    Y = d$Y, X = d$X, dName = d$dName,
    smp_data = d$smp, pop_data = d$pop,
    meanOnly = FALSE, MSE = "none", num.trees = 50, mtry = 2,
    aggregate_to = "subunit"
  )

  expect_equal(lme4::fixef(plain$MERFmodel$EffectModel),
               lme4::fixef(finer$MERFmodel$EffectModel), tolerance = 1e-8)
  expect_equal(plain$MERFmodel$RanEffSD, finer$MERFmodel$RanEffSD, tolerance = 1e-8)
  expect_equal(as.data.frame(plain$MERFmodel$VarianceCovariance),
               as.data.frame(finer$MERFmodel$VarianceCovariance), tolerance = 1e-8)
})

test_that("aggregate_to produces one indicator row per output group, coherent with the plain fit", {
  d <- tiny_saef_data_with_subunit()

  set.seed(2)
  finer <- SAEforest_model(
    Y = d$Y, X = d$X, dName = d$dName,
    smp_data = d$smp, pop_data = d$pop,
    meanOnly = FALSE, MSE = "none", num.trees = 50, mtry = 2,
    aggregate_to = "subunit"
  )

  expect_equal(nrow(finer$Indicators), length(unique(d$pop$subunit)))

  set.seed(2)
  plain <- SAEforest_model(
    Y = d$Y, X = d$X, dName = d$dName,
    smp_data = d$smp, pop_data = d$pop,
    meanOnly = FALSE, MSE = "none", num.trees = 50, mtry = 2
  )

  finer$Indicators$district <- sub("_[12]$", "", as.character(finer$Indicators$subunit))
  wts <- table(d$pop$subunit)
  finer$Indicators$w <- as.numeric(wts[as.character(finer$Indicators$subunit)])
  agg <- stats::aggregate(Mean * w ~ district, data = finer$Indicators, FUN = sum)
  wsum <- stats::aggregate(w ~ district, data = finer$Indicators, FUN = sum)
  agg_mean <- agg$`Mean * w` / wsum$w
  names(agg_mean) <- wsum$district

  plain_mean <- plain$Indicators$Mean
  names(plain_mean) <- as.character(plain$Indicators[[d$dName]])

  common <- intersect(names(agg_mean), names(plain_mean))
  expect_gt(length(common), 1)
  expect_equal(unname(agg_mean[common]), unname(plain_mean[common]), tolerance = 0.1)
})
