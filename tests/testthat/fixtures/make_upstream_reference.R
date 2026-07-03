# FROZEN baseline reference. Run ONCE, at the fork baseline, BEFORE any R/ change.
# Do not regenerate casually: this is the target the var.adjust=FALSE parity test
# must keep reproducing.
devtools::load_all(quiet = TRUE)
source(testthat::test_path("helper-fixtures.R"))
d <- tiny_saef_data()
set.seed(2026)
ref <- SAEforest_model(
  Y = d$Y, X = d$X, dName = d$dName,
  smp_data = d$smp, pop_data = d$pop,
  meanOnly = FALSE, MSE = "none",
  num.trees = 50, mtry = 3
)
saveRDS(ref$Indicators,
        testthat::test_path("fixtures", "upstream_nonLin_none.rds"))
