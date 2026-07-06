# Golden reference for the mean-path nonparametric MSE bootstrap, after fixing
# the bug where MSE_SAEforest_mean(...) in SAEforest_mean.R dropped the
# caller's `...` (num.trees, mtry, etc never reached the mean-path MSE
# bootstrap refits, which silently fell back to ranger's defaults). Run ONCE,
# after the fix, to freeze MSE_Estimates under the tiny fixture's forest
# config (num.trees = 50, mtry = 3). Guarded by test-mean-mse-dots.R.
devtools::load_all(quiet = TRUE)
source(testthat::test_path("helper-fixtures.R"))
d <- tiny_saef_data()
set.seed(2026)
res <- suppressMessages(SAEforest_model(
  Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp, pop_data = d$pop,
  meanOnly = TRUE, MSE = "nonparametric", B = 2, num.trees = 50, mtry = 3
))
saveRDS(res$MSE_Estimates,
        testthat::test_path("fixtures", "mean_mse_none.rds"))
